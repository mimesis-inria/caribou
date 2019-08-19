#ifdef CARIBOU_WITH_OPENMP
#include <omp.h>
#endif

#include <sofa/core/visual/VisualParams.h>
#include <sofa/core/ObjectFactory.h>
#include <sofa/simulation/Node.h>
#include <sofa/helper/AdvancedTimer.h>

#include <Caribou/Geometry/Tetrahedron.h>
#include <Caribou/Mechanics/Elasticity/Strain.h>

#include "TetrahedronElasticForce.h"

namespace SofaCaribou::GraphComponents::forcefield {
using namespace caribou::mechanics;

template<typename CanonicalTetrahedron>
TetrahedronElasticForce<CanonicalTetrahedron>::TetrahedronElasticForce()
: d_youngModulus(initData(&d_youngModulus,
    Real(1000), "youngModulus",
    "Young's modulus of the material",
    true /*displayed_in_GUI*/, false /*read_only_in_GUI*/))
, d_poissonRatio(initData(&d_poissonRatio,
    Real(0.3),  "poissonRatio",
    "Poisson's ratio of the material",
    true /*displayed_in_GUI*/, false /*read_only_in_GUI*/))
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

template<typename CanonicalTetrahedron>
void TetrahedronElasticForce<CanonicalTetrahedron>::init()
{
    Inherit::init();
    if (not d_topology_container.get()) {
        auto containers = this->getContext()->template getObjects<BaseMeshTopology>(BaseContext::Local);
        auto node = dynamic_cast<const sofa::simulation::Node *> (this->getContext());
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

template<typename CanonicalTetrahedron>
void TetrahedronElasticForce<CanonicalTetrahedron>::reinit()
{
    sofa::core::topology::BaseMeshTopology * topology = d_topology_container.get();
    MechanicalState<DataTypes> * state = this->mstate.get();

    if (!topology or !state)
        return;

    if (topology->getNbTetrahedra() == 0) {
        msg_warning() << "The topology container ('" << topology->getPathName() << "') does not contain any tetrahedron.";
        return;
    }

    const sofa::helper::ReadAccessor<Data<VecCoord>> X = state->readRestPositions();

    p_stiffness_matrices.resize(0);
    p_quadrature_nodes.resize(0);

    // Make sure every node of the tetrahedron have its coordinates inside the mechanical state vector
    for (std::size_t tetrahedron_id = 0; tetrahedron_id < topology->getNbTetrahedra(); ++tetrahedron_id) {
        const auto & node_indices = topology->getTetrahedron(tetrahedron_id);
        for (std::size_t j = 0; j < NumberOfNodes; ++j) {
            const auto & node_id = node_indices[j];
            if (node_id > X.size()-1) {
                msg_error() << "Some tetrahedrons have node indices outside of the state's position vector. Make sure "
                               "that the mechanical object '" << state->getPathName() << "' contains the position of "
                                                                                         "every tetrahedron nodes.";
                return;
            }
        }
    }

    p_initial_rotation.resize(topology->getNbTetrahedra(), Mat33::Identity());
    p_current_rotation.resize(topology->getNbTetrahedra(), Mat33::Identity());

    // Initialize the initial frame of each tetrahedron
    if (d_corotated.getValue()) {
        if (not d_linear_strain.getValue()) {
            msg_warning() << "The corotated method won't be computed since nonlinear strain is used.";
        } else {
            for (std::size_t tetrahedron_id = 0; tetrahedron_id < topology->getNbTetrahedra(); ++tetrahedron_id) {
                Tetrahedron tetra = tetrahedron(tetrahedron_id, X);
                p_initial_rotation[tetrahedron_id] = tetra.frame();
            }
        }
    }

    // Gather the integration points for each tetrahedron
    p_quadrature_nodes.resize(topology->getNbTetrahedra());
    for (std::size_t tetrahedron_id = 0; tetrahedron_id < topology->getNbTetrahedra(); ++tetrahedron_id) {
        auto   tetra = tetrahedron(tetrahedron_id, X);
        auto & quadrature_nodes = p_quadrature_nodes[tetrahedron_id];
        quadrature_nodes.resize(Tetrahedron::number_of_gauss_nodes);
        for (std::size_t gauss_node_id = 0; gauss_node_id < Tetrahedron::number_of_gauss_nodes; ++gauss_node_id) {
            const auto gauss_node  = MapVector<3>(Tetrahedron::gauss_nodes[gauss_node_id]);
            const auto &gauss_weight = Tetrahedron::gauss_weights[gauss_node_id];

            // Jacobian of the gauss node's transformation mapping from the elementary space to the world space
            const auto J = tetra.jacobian(gauss_node);
            const Mat33 Jinv = J.inverse();
            const auto detJ = J.determinant();

            // Derivatives of the shape functions at the gauss node with respect to global coordinates x,y and z
            const  Matrix<NumberOfNodes, 3> dN_dx = (Jinv.transpose() * Tetrahedron::dL(gauss_node).transpose()).transpose();
            quadrature_nodes[gauss_node_id] = {
                    gauss_weight,
                    detJ,
                    dN_dx,
                    Mat33::Identity()
            };
        }
    }

    // Initialize the stiffness matrix of every tetrahedrons
    p_stiffness_matrices.resize(topology->getNbTetrahedra());

    // Compute the initial tangent stiffness matrix
    compute_K();
}

template<typename CanonicalTetrahedron>
void TetrahedronElasticForce<CanonicalTetrahedron>::addForce(
        const MechanicalParams* mparams,
        Data<VecDeriv>& d_f,
        const Data<VecCoord>& d_x,
        const Data<VecDeriv>& d_v)
{
    SOFA_UNUSED(mparams);
    SOFA_UNUSED(d_v);

    static const auto I = Matrix<3,3>::Identity();

    auto topology = d_topology_container.get();
    MechanicalState<DataTypes> * state = this->mstate.get();

    if (!topology or !state)
        return;

    if (p_stiffness_matrices.size() != topology->getNbTetrahedra())
        return;

    sofa::helper::ReadAccessor<Data<VecCoord>> x = d_x;
    sofa::helper::ReadAccessor<Data<VecCoord>> x0 =  state->readRestPositions();
    sofa::helper::WriteAccessor<Data<VecDeriv>> f = d_f;

    const std::vector<Mat33> & initial_rotation = p_initial_rotation;
    std::vector<Mat33> & current_rotation = p_current_rotation;

    bool corotated = d_corotated.getValue();
    bool linear = d_linear_strain.getValue();
    recompute_compute_tangent_stiffness = (not linear) and mparams->implicit();
    if (linear) {
        // Small (linear) strain
        sofa::helper::AdvancedTimer::stepBegin("TetrahedronElasticForce::addForce");
        #pragma omp parallel for default(none) shared(topology, current_rotation, corotated, x, x0, f)
        for (std::size_t element_id = 0; element_id < topology->getNbTetrahedra(); ++element_id) {
            Tetrahedron tetra = tetrahedron(element_id, x);

            const Mat33 & R0 = initial_rotation[element_id];
            const Mat33 R0t = R0.transpose();

            Mat33 & R = current_rotation[element_id];


            // Extract the tetrahedron's frame
            if (corotated)
                R = tetra.frame();

            const Mat33 & Rt = R.transpose();

            // Gather the displacement vector
            Vector<NumberOfNodes*3> U;
            size_t i = 0;
            for (const auto &node_id : topology->getTetrahedron(element_id)) {
                const Vec3 r0 {x0[node_id][0], x0[node_id][1],  x0[node_id][2]};
                const Vec3 r  {x [node_id][0],  x [node_id][1], x [node_id][2]};

                const Vec3 u = Rt*r - R0t*r0;

                U[i++] = u[0];
                U[i++] = u[1];
                U[i++] = u[2];
            }

            // Compute the force vector
            const auto &K = p_stiffness_matrices[element_id];
            Vector<NumberOfNodes*3> F = K * U;

            // Write the forces into the output vector
            i = 0;
            for (const auto &node_id : topology->getTetrahedron(element_id)) {
                Vec3 force {F[i*3+0], F[i*3+1], F[i*3+2]};
                force = R*force;

                #pragma omp atomic
                f[node_id][0] -= force[0];

                #pragma omp atomic
                f[node_id][1] -= force[1];

                #pragma omp atomic
                f[node_id][2] -= force[2];
                ++i;
            }
        }
        sofa::helper::AdvancedTimer::stepEnd("TetrahedronElasticForce::addForce");
    } else {
        // Nonlinear Green-Lagrange strain
        sofa::helper::AdvancedTimer::stepBegin("TetrahedronElasticForce::addForce");
        const Real youngModulus = d_youngModulus.getValue();
        const Real poissonRatio = d_poissonRatio.getValue();

        const Real l = youngModulus * poissonRatio / ((1 + poissonRatio) * (1 - 2 * poissonRatio));
        const Real m = youngModulus / (2 * (1 + poissonRatio));

        #pragma omp parallel for default(none) shared(topology, current_rotation, corotated, x, x0, f, l, m)
        for (std::size_t element_id = 0; element_id < topology->getNbTetrahedra(); ++element_id) {
            const auto &tetra = topology->getTetrahedron(element_id);

            Matrix<NumberOfNodes, 3, Eigen::RowMajor> U;

            for (size_t i = 0; i < NumberOfNodes; ++i) {
                const auto u = x[tetra[i]] - x0[tetra[i]];
                U(i, 0) = u[0];
                U(i, 1) = u[1];
                U(i, 2) = u[2];
            }

            Matrix<NumberOfNodes, 3> forces;
            forces.fill(0);
            for (GaussNode &gauss_node : p_quadrature_nodes[element_id]) {

                // Jacobian of the gauss node's transformation mapping from the elementary space to the world space
                const auto detJ = gauss_node.jacobian_determinant;

                // Derivatives of the shape functions at the gauss node with respect to global coordinates x,y and z
                const auto dN_dx = gauss_node.dN_dx;

                // Gauss quadrature node weight
                const auto w = gauss_node.weight;

                // Deformation tensor at gauss node
                auto & F = gauss_node.F;
                F = elasticity::strain::F(dN_dx, U);

                // Strain tensor at gauss node
                const Mat33 C = F * F.transpose();
                const Mat33 E = 1/2. * (C - I);

                // Stress tensor at gauss node
                const Mat33 S = 2.*m*E + (l * E.trace() * I);

                // Elastic forces w.r.t the gauss node applied on each nodes
                for (size_t i = 0; i < NumberOfNodes; ++i) {
                    const Vec3 dx = dN_dx.row(i).transpose();
                    const Vec3 f_ = (detJ * w) * F.transpose() * (S * dx);
                    forces(i, 0) += f_[0];
                    forces(i, 1) += f_[1];
                    forces(i, 2) += f_[2];
                }
            }

            for (size_t i = 0; i < NumberOfNodes; ++i) {
                #pragma omp atomic
                f[tetra[i]][0] -= forces.row(i)[0];

                #pragma omp atomic
                f[tetra[i]][1] -= forces.row(i)[1];

                #pragma omp atomic
                f[tetra[i]][2] -= forces.row(i)[2];
            }
        }
        sofa::helper::AdvancedTimer::stepEnd("TetrahedronElasticForce::addForce");
    }
}

template<typename CanonicalTetrahedron>
void TetrahedronElasticForce<CanonicalTetrahedron>::addDForce(
        const MechanicalParams* mparams,
        Data<VecDeriv>& d_df,
        const Data<VecDeriv>& d_dx)
{
    auto * topology = d_topology_container.get();
    MechanicalState<DataTypes> * state = this->mstate.get();

    if (!topology or !state)
        return;

    if (p_stiffness_matrices.size() != topology->getNbTetrahedra())
        return;

    if (recompute_compute_tangent_stiffness)
        compute_K();

    auto kFactor = (Real)mparams->kFactorIncludingRayleighDamping(this->rayleighStiffness.getValue());
    sofa::helper::ReadAccessor<Data<VecDeriv>> dx = d_dx;
    sofa::helper::WriteAccessor<Data<VecDeriv>> df = d_df;
    std::vector<Mat33> & current_rotation = p_current_rotation;

    sofa::helper::AdvancedTimer::stepBegin("TetrahedronElasticForce::addDForce");
    #pragma omp parallel for default (none) shared (topology, current_rotation, dx, df, kFactor)
    for (std::size_t element_id = 0; element_id < topology->getNbTetrahedra(); ++element_id) {

        const Mat33 & R  = current_rotation[element_id];
        const Mat33 & Rt = R.transpose();

        // Gather the displacement vector
        Vector<NumberOfNodes*3> U;
        size_t i = 0;
        for (const auto & node_id : topology->getTetrahedron(element_id)) {
            const Vec3 v = {dx[node_id][0], dx[node_id][1], dx[node_id][2]};
            const Vec3 u = Rt*v;

            U[i++] = u[0];
            U[i++] = u[1];
            U[i++] = u[2];
        }

        // Compute the force vector
        const auto & K = p_stiffness_matrices[element_id];
        Vector<NumberOfNodes*3> F = K*U*kFactor;

        // Write the forces into the output vector
        i = 0;
        for (const auto & node_id : topology->getTetrahedron(element_id)) {
            Vec3 force {F[i*3+0], F[i*3+1], F[i*3+2]};
            force = R*force;

            #pragma omp atomic
            df[node_id][0] -= force[0];

            #pragma omp atomic
            df[node_id][1] -= force[1];

            #pragma omp atomic
            df[node_id][2] -= force[2];

            ++i;
        }
    }
    sofa::helper::AdvancedTimer::stepEnd("TetrahedronElasticForce::addDForce");
}

template<typename CanonicalTetrahedron>
void TetrahedronElasticForce<CanonicalTetrahedron>::addKToMatrix(
        sofa::defaulttype::BaseMatrix * matrix,
        SReal kFact,
        unsigned int & /*offset*/)
{
    auto * topology = d_topology_container.get();

    if (!topology)
        return;

    if (recompute_compute_tangent_stiffness)
        compute_K();

    std::vector<Mat33> & current_rotation = p_current_rotation;

    sofa::helper::AdvancedTimer::stepBegin("TetrahedronElasticForce::addKToMatrix");

    #pragma omp parallel for default (none) shared (topology, current_rotation, matrix, kFact)
    for (std::size_t element_id = 0; element_id < topology->getNbTetrahedra(); ++element_id) {
        const auto & node_indices = topology->getTetrahedron(element_id);
        const Mat33 & R  = current_rotation[element_id];
        const Mat33   Rt = R.transpose();

        const auto & K = p_stiffness_matrices[element_id];

        for (size_t i = 0; i < NumberOfNodes; ++i) {
            for (size_t j = 0; j < NumberOfNodes; ++j) {
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

                        #pragma omp critical
                        matrix->add(x, y, k(m,n));
                    }
                }
            }
        }
    }
    sofa::helper::AdvancedTimer::stepEnd("TetrahedronElasticForce::addKToMatrix");
}

template<typename CanonicalTetrahedron>
void TetrahedronElasticForce<CanonicalTetrahedron>::computeBBox(const sofa::core::ExecParams* params, bool onlyVisible)
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

template<typename CanonicalTetrahedron>
void TetrahedronElasticForce<CanonicalTetrahedron>::draw(const sofa::core::visual::VisualParams* vparams)
{
    auto *topology = d_topology_container.get();
    if (!topology)
        return;

    if (!vparams->displayFlags().getShowForceFields())
        return;

    vparams->drawTool()->saveLastState();

    if (vparams->displayFlags().getShowWireFrame())
        vparams->drawTool()->setPolygonMode(0, true);

    vparams->drawTool()->disableLighting();

    const VecCoord &x = this->mstate->read(sofa::core::ConstVecCoordId::position())->getValue();

    std::vector< sofa::defaulttype::Vec<3, double> > points[4];
    for(Topology::TetrahedronID i = 0 ; i<topology->getNbTetrahedra();++i)
    {
        const auto t=topology->getTetra(i);

        const auto & a = t[0];
        const auto & b = t[1];
        const auto & c = t[2];
        const auto & d = t[3];
        Coord center = (x[a]+x[b]+x[c]+x[d])*0.125;
        Coord pa = (x[a]+center)*(Real)0.666667;
        Coord pb = (x[b]+center)*(Real)0.666667;
        Coord pc = (x[c]+center)*(Real)0.666667;
        Coord pd = (x[d]+center)*(Real)0.666667;

        points[0].push_back(pa);
        points[0].push_back(pb);
        points[0].push_back(pc);

        points[1].push_back(pb);
        points[1].push_back(pc);
        points[1].push_back(pd);

        points[2].push_back(pc);
        points[2].push_back(pd);
        points[2].push_back(pa);

        points[3].push_back(pd);
        points[3].push_back(pa);
        points[3].push_back(pb);
    }

    sofa::defaulttype::Vec<4, float> face_colors[4] = {
            {1.0, 0.0, 0.0, 1.0},
            {1.0, 0.0, 0.5, 1.0},
            {1.0, 1.0, 0.0, 1.0},
            {1.0, 0.5, 1.0, 1.0}
    };

    vparams->drawTool()->drawTriangles(points[0], face_colors[0]);
    vparams->drawTool()->drawTriangles(points[1], face_colors[1]);
    vparams->drawTool()->drawTriangles(points[2], face_colors[2]);
    vparams->drawTool()->drawTriangles(points[3], face_colors[3]);

    if (vparams->displayFlags().getShowWireFrame())
        vparams->drawTool()->setPolygonMode(0,false);

    vparams->drawTool()->restoreLastState();
}

template<typename CanonicalTetrahedron>
void TetrahedronElasticForce<CanonicalTetrahedron>::compute_K()
{
    static const auto I = Matrix<3,3>::Identity();
    auto topology = d_topology_container.get();

    if (!topology)
        return;

    if (p_stiffness_matrices.size() != topology->getNbTetrahedra())
        return;

    const Real youngModulus = d_youngModulus.getValue();
    const Real poissonRatio = d_poissonRatio.getValue();

    const Real l = youngModulus * poissonRatio / ((1 + poissonRatio) * (1 - 2 * poissonRatio));
    const Real m = youngModulus / (2 * (1 + poissonRatio));

    sofa::helper::AdvancedTimer::stepBegin("TetrahedronElasticForce::compute_k");

    #pragma omp parallel for default (none) shared (topology, l, m)
    for (std::size_t element_id = 0; element_id < topology->getNbTetrahedra(); ++element_id) {
        auto & K = p_stiffness_matrices[element_id];
        K.fill(0.);

        for (GaussNode &gauss_node : p_quadrature_nodes[element_id]) {
            // Jacobian of the gauss node's transformation mapping from the elementary space to the world space
            const auto detJ = gauss_node.jacobian_determinant;

            // Derivatives of the shape functions at the gauss node with respect to global coordinates x,y and z
            const auto dN_dx = gauss_node.dN_dx;

            // Gauss quadrature node weight
            const auto w = gauss_node.weight;

            // Deformation tensor at gauss node
            auto & F = gauss_node.F;

            // Strain tensor at gauss node
            const auto C = F * F.transpose();
            const auto E = 1/2. * (C - I);

            // Stress tensor at gauss node
            const auto S = 2.*m*E + (l * E.trace() * I);

            // Computation of the tangent-stiffness matrix
            for (std::size_t i = 0; i < NumberOfNodes; ++i) {
                // Derivatives of the ith shape function at the gauss node with respect to global coordinates x,y and z
                const auto dxi = dN_dx.row(i).transpose();

                for (std::size_t j = 0; j < NumberOfNodes; ++j) {
                    // Derivatives of the jth shape function at the gauss node with respect to global coordinates x,y and z
                    const auto dxj = dN_dx.row(j).transpose();

                    // Derivative of the force applied on node j w.r.t the u component of the ith nodal's displacement
                    const auto dFu = dxi * I.row(0); // Deformation tensor derivative with respect to u_i
                    const auto dCu = dFu*F.transpose() + F*dFu.transpose();
                    const auto dEu = 1/2. * dCu;
                    const auto dSu = 2. * m * dEu + (l*dEu.trace() * I);
                    const Vec3  Ku  = (detJ*w) * (dFu.transpose()*S + F.transpose()*dSu) * dxj;

                    // Derivative of the force applied on node j w.r.t the v component of the ith nodal's displacement
                    const auto dFv = dxi * I.row(1); // Deformation tensor derivative with respect to u_i
                    const auto dCv = dFv*F.transpose() + F*dFv.transpose();
                    const auto dEv = 1/2. * dCv;
                    const auto dSv = 2. * m * dEv + (l*dEv.trace() * I);
                    const Vec3  Kv  = (detJ*w) * (dFv.transpose()*S + F.transpose()*dSv) * dxj;

                    // Derivative of the force applied on node j w.r.t the w component of the ith nodal's displacement
                    const auto dFw = dxi * I.row(2); // Deformation tensor derivative with respect to u_i
                    const auto dCw = dFw*F.transpose() + F*dFw.transpose();
                    const auto dEw = 1/2. * dCw;
                    const auto dSw = 2. * m * dEw + (l*dEw.trace() * I);
                    const Vec3  Kw  = (detJ*w) * (dFw.transpose()*S + F.transpose()*dSw) * dxj;

                    Mat33 Kji;
                    Kji << Ku, Kv, Kw; // Kji * ui = fj (the force applied on node j on the displacement of the node i)

                    for (size_t ii = 0; ii < 3; ++ii) {
                        for (size_t jj = 0; jj < 3; ++jj) {
                            K(j*3 +ii, i*3 + jj) += Kji(ii,jj);
                        }
                    }
                }
            }
        }
    }
    recompute_compute_tangent_stiffness = false;
    K_is_up_to_date = false;
    eigenvalues_are_up_to_date = false;
    sofa::helper::AdvancedTimer::stepEnd("TetrahedronElasticForce::compute_k");
}

static int TetrahedronElasticForceClass = RegisterObject("Caribou tetrahedron FEM Forcefield")
    .add< TetrahedronElasticForce<caribou::geometry::interpolation::Tetrahedron4>>(true)
;

} // namespace SofaCaribou::GraphComponents::forcefield