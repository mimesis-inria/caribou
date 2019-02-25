#include <numeric>

#include <sofa/core/visual/VisualParams.h>
#include <sofa/core/ObjectFactory.h>
#include <sofa/simulation/Node.h>
#include <sofa/helper/AdvancedTimer.h>

#include <SofaCaribou/Traits.h>
#include <Caribou/Geometry/Hexahedron.h>
#include <Caribou/Mechanics/Elasticity/Strain.h>

#include "HexahedronElasticForce.h"


namespace SofaCaribou::GraphComponents::forcefield {

using namespace sofa::core::topology;
using namespace caribou::geometry;
using namespace caribou::algebra;
using namespace caribou::mechanics;

template<class DataTypes>
HexahedronElasticForce<DataTypes>::HexahedronElasticForce()
: d_youngModulus(initData(&d_youngModulus,
        Real(1000), "youngModulus", "Young's modulus of the material", true /*displayed_in_GUI*/, false /*read_only_in_GUI*/))
, d_poissonRatio(initData(&d_poissonRatio,
        Real(0.3),  "poissonRatio", "Poisson's ratio of the material", true /*displayed_in_GUI*/, false /*read_only_in_GUI*/))
, d_linear_strain(initData(&d_linear_strain,
        bool(true), "linearStrain",
        "True if the small (linear) strain tensor is used, otherwise the nonlinear Green-Lagrange strain is used.",
        true /*displayed_in_GUI*/, false /*read_only_in_GUI*/))
, d_corotated(initData(&d_corotated,
        bool(true), "corotated",
        "Whether or not to use corotated elements for the strain computation.",
        true /*displayed_in_GUI*/, false /*read_only_in_GUI*/))
, d_topology_container(initLink(
        "topology_container", "Topology that contains the elements on which this force will be computed."))
{
}

template<class DataTypes>
void HexahedronElasticForce<DataTypes>::init()
{
    Inherit::init();
    if (not d_topology_container.get()) {
        auto containers = this->getContext()->template getObjects<BaseMeshTopology>(BaseContext::Local);
        auto node = static_cast<const sofa::simulation::Node *> (this->getContext());
        if (containers.empty()) {
            msg_error() << "No topology were found in the context node '" << node->getPathName() << "'.";
        } else if (containers.size() > 1) {
            msg_error() <<
            "Multiple topology were found in the node '" << node->getPathName() << "'." <<
            " Please specify which one contains the elements on which this force field will be computed " <<
            "by explicitly setting the container's path in the 'topology_container' parameter.";
        } else {
            d_topology_container.set(containers[0]);
            msg_info() << "Automatically found the topology '" << d_topology_container.get()->getPathName() << "'.";
        }
    }

    if (!this->mstate.get())
        msg_error() << "No mechanical state set. Please add a mechanical object in the current node or specify a path "
                       "to one in the 'mstate' parameter.";


    reinit();
}

template<class DataTypes>
void HexahedronElasticForce<DataTypes>::reinit()
{
    sofa::core::topology::BaseMeshTopology * topology = d_topology_container.get();
    MechanicalState<DataTypes> * state = this->mstate.get();

    if (!topology or !state)
        return;

    if (topology->getNbHexahedra() == 0) {
        msg_warning() << "The topology container ('" << topology->getPathName() << "') does not contain any hexahedron.";
        return;
    }

    const sofa::helper::ReadAccessor<Data<VecCoord>> X = state->readRestPositions();

    p_stiffness_matrices.resize(0);
    p_quatrature_nodes.resize(0);

    // Make sure every node of the hexahedrons have its coordinates inside the mechanical state vector
    for (std::size_t hexa_id = 0; hexa_id < topology->getNbHexahedra(); ++hexa_id) {
        const auto & node_indices = topology->getHexahedron(hexa_id);
        for (std::size_t j = 0; j < 8; ++j) {
            const auto & node_id = node_indices[j];
            if (node_id > X.size()-1) {
                msg_error() << "Some hexahedrons have node indices outside of the state's position vector. Make sure "
                               "that the mechanical object '" << state->getPathName() << "' contains the position of "
                                                                                         "every hexahedron nodes.";
                return;
            }
        }
    }

    p_initial_rotation.resize(topology->getNbHexahedra(), Mat33::Identity());
    p_current_rotation.resize(topology->getNbHexahedra(), Mat33::Identity());

    // Initialize the initial frame of each hexahedron
    if (d_corotated.getValue()) {
        if (not d_linear_strain.getValue()) {
            msg_warning() << "The corotated method won't be computed since nonlinear strain is used.";
        } else {
            for (std::size_t hexa_id = 0; hexa_id < topology->getNbHexahedra(); ++hexa_id) {
                Hexahedron hexa = make_hexa(hexa_id, X);
                p_initial_rotation[hexa_id] = hexa.extract_frame_by_cross_products();
            }
        }
    }

    if (not d_linear_strain.getValue()) {
        // For nonlinear Green-Lagrange strain
        // Initialize the gauss quadrature points for every hexahedrons
        p_quatrature_nodes.resize(topology->getNbHexahedra());
        for (std::size_t hexa_id = 0; hexa_id < topology->getNbHexahedra(); ++hexa_id) {
            Hexahedron hexa = make_hexa(hexa_id, X);

            for (std::size_t gauss_node_id = 0; gauss_node_id < Hexahedron::gauss_nodes.size(); ++gauss_node_id) {
                const auto &gauss_node = Hexahedron::gauss_nodes[gauss_node_id];

                // Local coordinates of the gauss node
                const auto &u = gauss_node[0];
                const auto &v = gauss_node[1];
                const auto &w = gauss_node[2];

                // Jacobian of the gauss node's transformation mapping from the elementary space to the world space
                const auto J = hexa.jacobian(u, v, w);
                const auto Jinv = J.inverted();
                const auto detJ = J.determinant();

                // Derivatives of the shape functions at the gauss node with respect to global coordinates x,y and z
                const auto dN_dx = (Jinv.T() * interpolation::Hexahedron8::dN({u, v, w}).T()).T();

                p_quatrature_nodes[hexa_id][gauss_node_id] = {
                        Hexahedron::gauss_weights[gauss_node_id],
                        detJ,
                        dN_dx
                };
            }

        }
    }

    // Initialize the stiffness matrix of every hexahedrons
    p_stiffness_matrices.resize(topology->getNbHexahedra());

    const Real youngModulus = d_youngModulus.getValue();
    const Real poissonRatio = d_poissonRatio.getValue();

    const Real l = youngModulus*poissonRatio / ((1 + poissonRatio)*(1 - 2*poissonRatio));
    const Real m = youngModulus / (2 * (1 + poissonRatio));
    const Real a = l + 2*m;
    Matrix<6,6,Real> C;
    C(0,0) = a; C(0,1) = l; C(0,2) = l; C(0,3) = 0; C(0,4) = 0; C(0,5) = 0;
    C(1,0) = l; C(1,1) = a; C(1,2) = l; C(1,3) = 0; C(1,4) = 0; C(1,5) = 0;
    C(2,0) = l; C(2,1) = l; C(2,2) = a; C(2,3) = 0; C(2,4) = 0; C(2,5) = 0;
    C(3,0) = 0; C(3,1) = 0; C(3,2) = 0; C(3,3) = m; C(3,4) = 0; C(3,5) = 0;
    C(4,0) = 0; C(4,1) = 0; C(4,2) = 0; C(4,3) = 0; C(4,4) = m; C(4,5) = 0;
    C(5,0) = 0; C(5,1) = 0; C(5,2) = 0; C(5,3) = 0; C(5,4) = 0; C(5,5) = m;

    for (std::size_t hexa_id = 0; hexa_id < topology->getNbHexahedra(); ++hexa_id) {
        Hexahedron hexa = make_hexa(hexa_id, X);
        p_stiffness_matrices[hexa_id] = hexa.gauss_quadrature([&C](const auto & hexa, const auto & u, const auto & v, const auto & w) {
            const auto B = elasticity::strain::B(hexa, u, v, w);
            return B.T() * C * B;
        });
    }
}

template<class DataTypes>
void HexahedronElasticForce<DataTypes>::addForce(
        const MechanicalParams * mparams,
        Data<VecDeriv> & d_f,
        const Data<VecCoord> & d_x,
        const Data<VecDeriv> & d_v)
{
    SOFA_UNUSED(mparams);
    SOFA_UNUSED(d_v);

    static const auto I = Matrix<3,3,Real>::Identity();

    auto topology = d_topology_container.get();
    MechanicalState<DataTypes> * state = this->mstate.get();

    if (!topology or !state)
        return;

    if (p_stiffness_matrices.size() != topology->getNbHexahedra())
        return;

    sofa::helper::ReadAccessor<Data<VecCoord>> x = d_x;
    sofa::helper::ReadAccessor<Data<VecCoord>> x0 =  state->readRestPositions();
    sofa::helper::WriteAccessor<Data<VecDeriv>> f = d_f;

    const std::vector<Mat33> & initial_rotation = p_initial_rotation;
    std::vector<Mat33> & current_rotation = p_current_rotation;

    bool corotated = d_corotated.getValue();
    if (d_linear_strain.getValue()) {
        // Small (linear) strain
        sofa::helper::AdvancedTimer::stepBegin("HexahedronElasticForce::addForce");
        for (std::size_t hexa_id = 0; hexa_id < topology->getNbHexahedra(); ++hexa_id) {
            Hexahedron hexa = make_hexa(hexa_id, x);

            const Mat33 & R0 = initial_rotation[hexa_id];
            const Mat33 R0t = R0.T();

            Mat33 & R = current_rotation[hexa_id];


            // Extract the hexahedron's frame
            if (corotated)
                R = hexa.extract_frame_by_cross_products();

            const Mat33 & Rt = R.T();

            // Gather the displacement vector
            Vec24 U;
            size_t i = 0;
            for (const auto &node_id : topology->getHexahedron(hexa_id)) {
                const Vec3 r0 {x0[node_id][0], x0[node_id][1],  x0[node_id][2]};
                const Vec3 r  {x [node_id][0],  x [node_id][1], x [node_id][2]};

                const Vec3 u = Rt*r - R0t*r0;

                U[i++] = u[0];
                U[i++] = u[1];
                U[i++] = u[2];
            }

            // Compute the force vector
            const auto &K = p_stiffness_matrices[hexa_id];
            Vec24 F = K * U;

            // Write the forces into the output vector
            i = 0;
            for (const auto &node_id : topology->getHexahedron(hexa_id)) {
                Vec3 force {F[i*3+0], F[i*3+1], F[i*3+2]};
                force = R*force;

                f[node_id][0] -= force[0];
                f[node_id][1] -= force[1];
                f[node_id][2] -= force[2];
                ++i;
            }
        }
        sofa::helper::AdvancedTimer::stepEnd("HexahedronElasticForce::addForce");
    } else {
        // Nonlinear Green-Lagrange strain
        sofa::helper::AdvancedTimer::stepBegin("HexahedronElasticForce::addForce");
        bool compute_tangent_stiffness = mparams->implicit();
        const Real youngModulus = d_youngModulus.getValue();
        const Real poissonRatio = d_poissonRatio.getValue();

        const Real l = youngModulus * poissonRatio / ((1 + poissonRatio) * (1 - 2 * poissonRatio));
        const Real m = youngModulus / (2 * (1 + poissonRatio));

        if (compute_tangent_stiffness) {
            for (Mat2424 & K : p_stiffness_matrices)
                K.fill(0.);
        }

        for (std::size_t hexa_id = 0; hexa_id < topology->getNbHexahedra(); ++hexa_id) {
            const auto &hexa = topology->getHexahedron(hexa_id);

            Matrix<8, 3> U;

            for (size_t i = 0; i < 8; ++i) {
                const auto u = x[hexa[i]] - x0[hexa[i]];
                U(i, 0) = u[0];
                U(i, 1) = u[1];
                U(i, 2) = u[2];
            }

            Matrix<8, 3, Real> forces;
            forces.fill(0);
            Mat2424 & K = p_stiffness_matrices[hexa_id];
            for (const GaussNode &gauss_node : p_quatrature_nodes[hexa_id]) {

                // Jacobian of the gauss node's transformation mapping from the elementary space to the world space
                const auto detJ = gauss_node.jacobian_determinant;

                // Derivatives of the shape functions at the gauss node with respect to global coordinates x,y and z
                const auto dN_dx = gauss_node.dN_dx;

                // Gauss quadrature node weight
                const auto w = gauss_node.weight;

                // Deformation tensor at gauss node
                const auto F = elasticity::strain::F(dN_dx, U);

                // Strain tensor at gauss node
                const auto C = F * F.T();
                const auto E = 1/2. * (C - I);

                // Stress tensor at gauss node
                const auto S = 2.*m*E + (l * tr(E) * I);

                // Elastic forces w.r.t the gauss node applied on each nodes
                for (size_t i = 0; i < 8; ++i) {
                    const auto dx = dN_dx.row(i).T();
                    const auto f_ = (detJ * w) * F.T() * (S * dx);
                    forces(i, 0) += f_[0];
                    forces(i, 1) += f_[1];
                    forces(i, 2) += f_[2];
                }

                // Computation of the tangent-stiffness matrix
                if (compute_tangent_stiffness) {
                    for (std::size_t i = 0; i < 8; ++i) {
                        for (std::size_t j = 0; j < 8; ++j) {

                            // Derivatives of the ith shape function at the gauss node with respect to global coordinates x,y and z
                            const auto dxi = dN_dx.row(i).T();

                            // Derivatives of the jth shape function at the gauss node with respect to global coordinates x,y and z
                            const auto dxj = dN_dx.row(j).T();

                            // Derivative of the force applied on node j w.r.t the u component of the ith nodal's displacement
                            const Mat33 dFu = dxi * I.row(0); // Deformation tensor derivative with respect to u_i
                            const Mat33 dCu = dFu*F.T() + F*dFu.T();
                            const Mat33 dEu = 1/2. * dCu;
                            const Mat33 dSu = 2. * m * dEu + (l*tr(dEu) * I);
                            const Vec3  Ku  = (detJ*w) * (dFu.T()*S + F.T()*dSu) * dxj;

                            // Derivative of the force applied on node j w.r.t the v component of the ith nodal's displacement
                            const Mat33 dFv = dxi * I.row(1); // Deformation tensor derivative with respect to u_i
                            const Mat33 dCv = dFv*F.T() + F*dFv.T();
                            const Mat33 dEv = 1/2. * dCv;
                            const Mat33 dSv = 2. * m * dEv + (l*tr(dEv) * I);
                            const Vec3  Kv  = (detJ*w) * (dFv.T()*S + F.T()*dSv) * dxj;

                            // Derivative of the force applied on node j w.r.t the w component of the ith nodal's displacement
                            const Mat33 dFw = dxi * I.row(2); // Deformation tensor derivative with respect to u_i
                            const Mat33 dCw = dFw*F.T() + F*dFw.T();
                            const Mat33 dEw = 1/2. * dCw;
                            const Mat33 dSw = 2. * m * dEw + (l*tr(dEw) * I);
                            const Vec3  Kw  = (detJ*w) * (dFw.T()*S + F.T()*dSw) * dxj;

                            const Mat33 Kji (Ku, Kv, Kw); // Kji * ui = fj (the force applied on node j on the displacement of the node i)

                            for (size_t ii = 0; ii < 3; ++ii) {
                                for (size_t jj = 0; jj < 3; ++jj) {
                                    K(j*3 +ii, i*3 + jj) += Kji(ii,jj);
                                }
                            }
                        }
                    }
                }
            }

            for (size_t i = 0; i < 8; ++i) {
                f[hexa[i]][0] -= forces.row(i)[0];
                f[hexa[i]][1] -= forces.row(i)[1];
                f[hexa[i]][2] -= forces.row(i)[2];
            }

        }
        sofa::helper::AdvancedTimer::stepEnd("HexahedronElasticForce::addForce");
    }
}

template<class DataTypes>
void HexahedronElasticForce<DataTypes>::addDForce(
        const MechanicalParams* mparams,
        Data<VecDeriv>& d_df,
        const Data<VecDeriv>& d_dx)
{
    auto * topology = d_topology_container.get();
    MechanicalState<DataTypes> * state = this->mstate.get();

    if (!topology or !state)
        return;

    if (p_stiffness_matrices.size() != topology->getNbHexahedra())
        return;

    auto kFactor = (Real)mparams->kFactorIncludingRayleighDamping(this->rayleighStiffness.getValue());
    sofa::helper::ReadAccessor<Data<VecDeriv>> dx = d_dx;
    sofa::helper::WriteAccessor<Data<VecDeriv>> df = d_df;
    std::vector<Mat33> & current_rotation = p_current_rotation;

    sofa::helper::AdvancedTimer::stepBegin("HexahedronElasticForce::addDForce");
    for (std::size_t hexa_id = 0; hexa_id < topology->getNbHexahedra(); ++hexa_id) {

        const Mat33 & R  = current_rotation[hexa_id];
        const Mat33 & Rt = R.T();

        // Gather the displacement vector
        Vec24 U;
        size_t i = 0;
        for (const auto & node_id : topology->getHexahedron(hexa_id)) {
            const Vec3 v = {dx[node_id][0], dx[node_id][1], dx[node_id][2]};
            const Vec3 u = Rt*v;

            U[i++] = u[0];
            U[i++] = u[1];
            U[i++] = u[2];
        }

        // Compute the force vector
        const auto & K = p_stiffness_matrices[hexa_id];
        Vec24 F = K*U*kFactor;

        // Write the forces into the output vector
        i = 0;
        for (const auto & node_id : topology->getHexahedron(hexa_id)) {
            Vec3 force {F[i*3+0], F[i*3+1], F[i*3+2]};
            force = R*force;

            df[node_id][0] -= force[0];
            df[node_id][1] -= force[1];
            df[node_id][2] -= force[2];
            ++i;
        }
    }
    sofa::helper::AdvancedTimer::stepEnd("HexahedronElasticForce::addDForce");
}

template<class DataTypes>
void HexahedronElasticForce<DataTypes>::addKToMatrix(sofa::defaulttype::BaseMatrix * matrix, SReal kFact, unsigned int & /*offset*/)
{
    auto * topology = d_topology_container.get();

    if (!topology)
        return;

    std::vector<Mat33> & current_rotation = p_current_rotation;

    sofa::helper::AdvancedTimer::stepBegin("HexahedronElasticForce::addKToMatrix");
    for (std::size_t hexa_id = 0; hexa_id < topology->getNbHexahedra(); ++hexa_id) {
        const auto & node_indices = topology->getHexahedron(hexa_id);
        const Mat33 & R  = current_rotation[hexa_id];
        const Mat33   Rt = R.T();

        const auto & K = p_stiffness_matrices[hexa_id];

        for (size_t i = 0; i < 8; ++i) {
            for (size_t j = 0; j < 8; ++j) {
                Mat33 k;

                for (unsigned char m = 0; m < 3; ++m) {
                    for (unsigned char n = 0; n < 3; ++n) {
                        k(m,n) = K(i*3+m, j*3+n);
                    }
                }

                k = -1. * R*k*Rt*kFact;

                for (unsigned char m = 0; m < 3; ++m) {
                    for (unsigned char n = 0; n < 3; ++n) {
                        const auto x = node_indices[i]*3+m;
                        const auto y = node_indices[j]*3+n;

                        matrix->add(x, y, k(m,n));
                    }
                }
            }
        }
    }
    sofa::helper::AdvancedTimer::stepEnd("HexahedronElasticForce::addKToMatrix");
}

template<class DataTypes>
void HexahedronElasticForce<DataTypes>::computeBBox(const sofa::core::ExecParams* params, bool onlyVisible)
{
    if( !onlyVisible ) return;

    sofa::helper::ReadAccessor<Data<VecCoord>> x = this->mstate->read(sofa::core::VecCoordId::position());

    static const Real max_real = std::numeric_limits<Real>::max();
    static const Real min_real = std::numeric_limits<Real>::lowest();
    Real maxBBox[3] = {min_real,min_real,min_real};
    Real minBBox[3] = {max_real,max_real,max_real};
    for (size_t i=0; i<x.size(); i++)
    {
        for (int c=0; c<3; c++)
        {
            if (x[i][c] > maxBBox[c]) maxBBox[c] = (Real)x[i][c];
            else if (x[i][c] < minBBox[c]) minBBox[c] = (Real)x[i][c];
        }
    }

    this->f_bbox.setValue(params,sofa::defaulttype::TBoundingBox<Real>(minBBox,maxBBox));
}

template<class DataTypes>
void HexahedronElasticForce<DataTypes>::draw(const sofa::core::visual::VisualParams* vparams)
{
    auto * topology = d_topology_container.get();
    if (!topology)
        return;

    if (!vparams->displayFlags().getShowForceFields())
        return;

    if (vparams->displayFlags().getShowWireFrame())
        vparams->drawTool()->setPolygonMode(0,true);

    vparams->drawTool()->disableLighting();

    const VecCoord& x = this->mstate->read(sofa::core::ConstVecCoordId::position())->getValue();


    for (auto node_indices : topology->getHexahedra()) {
        std::vector< sofa::defaulttype::Vector3 > points[6];

        auto a = node_indices[0];
        auto b = node_indices[1];
        auto d = node_indices[3];
        auto c = node_indices[2];
        auto e = node_indices[4];
        auto f = node_indices[5];
        auto h = node_indices[7];
        auto g = node_indices[6];


        Coord center = (x[a]+x[b]+x[c]+x[d]+x[e]+x[g]+x[f]+x[h])*0.125;
        Real percentage = 0.15;
        Coord pa = x[a]-(x[a]-center)*percentage;
        Coord pb = x[b]-(x[b]-center)*percentage;
        Coord pc = x[c]-(x[c]-center)*percentage;
        Coord pd = x[d]-(x[d]-center)*percentage;
        Coord pe = x[e]-(x[e]-center)*percentage;
        Coord pf = x[f]-(x[f]-center)*percentage;
        Coord pg = x[g]-(x[g]-center)*percentage;
        Coord ph = x[h]-(x[h]-center)*percentage;



        points[0].push_back(pa);
        points[0].push_back(pb);
        points[0].push_back(pc);
        points[0].push_back(pa);
        points[0].push_back(pc);
        points[0].push_back(pd);

        points[1].push_back(pe);
        points[1].push_back(pf);
        points[1].push_back(pg);
        points[1].push_back(pe);
        points[1].push_back(pg);
        points[1].push_back(ph);

        points[2].push_back(pc);
        points[2].push_back(pd);
        points[2].push_back(ph);
        points[2].push_back(pc);
        points[2].push_back(ph);
        points[2].push_back(pg);

        points[3].push_back(pa);
        points[3].push_back(pb);
        points[3].push_back(pf);
        points[3].push_back(pa);
        points[3].push_back(pf);
        points[3].push_back(pe);

        points[4].push_back(pa);
        points[4].push_back(pd);
        points[4].push_back(ph);
        points[4].push_back(pa);
        points[4].push_back(ph);
        points[4].push_back(pe);

        points[5].push_back(pb);
        points[5].push_back(pc);
        points[5].push_back(pg);
        points[5].push_back(pb);
        points[5].push_back(pg);
        points[5].push_back(pf);


        vparams->drawTool()->setLightingEnabled(false);
        vparams->drawTool()->drawTriangles(points[0], sofa::defaulttype::Vec<4,float>(0.7f,0.7f,0.1f,1.0f));
        vparams->drawTool()->drawTriangles(points[1], sofa::defaulttype::Vec<4,float>(0.7f,0.0f,0.0f,1.0f));
        vparams->drawTool()->drawTriangles(points[2], sofa::defaulttype::Vec<4,float>(0.0f,0.7f,0.0f,1.0f));
        vparams->drawTool()->drawTriangles(points[3], sofa::defaulttype::Vec<4,float>(0.0f,0.0f,0.7f,1.0f));
        vparams->drawTool()->drawTriangles(points[4], sofa::defaulttype::Vec<4,float>(0.1f,0.7f,0.7f,1.0f));
        vparams->drawTool()->drawTriangles(points[5], sofa::defaulttype::Vec<4,float>(0.7f,0.1f,0.7f,1.0f));

    }
}

} // namespace SofaCaribou::GraphComponents::forcefield

SOFA_DECL_CLASS(HexahedronElasticForce)
static int HexahedronElasticForceClass = sofa::core::RegisterObject("Caribou Hexahedron FEM Forcefield")
                                        .add< SofaCaribou::GraphComponents::forcefield::HexahedronElasticForce<sofa::defaulttype::Vec3dTypes> >(true)
;