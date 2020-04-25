#include <numeric>
#include <queue>
#include <array>

#include <sofa/core/visual/VisualParams.h>
#include <sofa/core/ObjectFactory.h>
#include <sofa/simulation/Node.h>
#include <sofa/helper/AdvancedTimer.h>

#include "HexahedronElasticForce.h"

#include <Eigen/Sparse>

#include <Caribou/Geometry/Hexahedron.h>
#include <Caribou/Mechanics/Elasticity/Strain.h>

#if !EIGEN_VERSION_AT_LEAST(3,3,0)
namespace Eigen {
using Index = EIGEN_DEFAULT_DENSE_INDEX_TYPE;
}
#endif


namespace SofaCaribou::forcefield {

using namespace sofa::core::topology;
using namespace caribou::geometry;
using namespace caribou::mechanics;

using sofa::defaulttype::Vec3Types;

HexahedronElasticForce::HexahedronElasticForce()
: d_youngModulus(initData(&d_youngModulus,
        Real(1000), "youngModulus",
        "Young's modulus of the material",
        true /*displayed_in_GUI*/, false /*read_only_in_GUI*/))
, d_poissonRatio(initData(&d_poissonRatio,
        Real(0.3),  "poissonRatio",
        "Poisson's ratio of the material",
        true /*displayed_in_GUI*/, false /*read_only_in_GUI*/))
, d_corotated(initData(&d_corotated,
        bool(true), "corotated",
        "Whether or not to use corotated elements for the strain computation.",
        true /*displayed_in_GUI*/, false /*read_only_in_GUI*/))
, d_integration_method(initData(&d_integration_method,
        "integration_method",
        R"(
                Integration method used to integrate the stiffness matrix.

                Methods are:
                  Regular:       Regular 8 points gauss integration (default).
                  OnePointGauss: One gauss point integration at the center of the hexahedron
                )",
        true /*displayed_in_GUI*/, false /*read_only_in_GUI*/))
, d_topology_container(initLink(
        "topology_container", "Topology that contains the elements on which this force will be computed."))
{
    d_integration_method.setValue(sofa::helper::OptionsGroup(std::vector<std::string> {
        "Regular", "OnePointGauss"
    }));

    sofa::helper::WriteAccessor<Data< sofa::helper::OptionsGroup >> integration_method = d_integration_method;
    integration_method->setSelectedItem(static_cast<unsigned int>(0));
}

void HexahedronElasticForce::init()
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
            "by explicitly setting the container's path in the  '" << d_topology_container.getName() << "'  parameter.";
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

void HexahedronElasticForce::reinit()
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

    // Gather the integration points for each hexahedron
    std::vector<std::vector<GaussNode>> quadrature_nodes(topology->getNbHexahedra());
    p_quadrature_nodes.resize(topology->getNbHexahedra());

    for (std::size_t hexa_id = 0; hexa_id < topology->getNbHexahedra(); ++hexa_id) {
        auto   hexa = hexahedron(hexa_id, X);

        // List of pair (Xi, w) where Xi is the local coordinates vector of the gauss point and w is its weight.
        std::vector<std::pair<Vec3, Real>> gauss_points;

        if (integration_method() == IntegrationMethod::OnePointGauss) {
            gauss_points.emplace_back(Vector<3>(0, 0, 0), 8);
        } else {
            for (const auto & gauss_node : hexa.gauss_nodes()) {
                const auto & position   = gauss_node.position;
                const auto & weight = gauss_node.weight;

                gauss_points.emplace_back(position, weight);
            }
        }

        // 2. At this point, we have all the gauss points of the hexa and their corrected weight.
        //    We now compute constant values required by the simulation for each of them
        auto &hexa_quadrature_nodes = p_quadrature_nodes[hexa_id];
        for (const auto &gauss_point : gauss_points) {
            // Jacobian of the gauss node's transformation mapping from the elementary space to the world space
            const auto J = hexa.jacobian(gauss_point.first);
            const Mat33 Jinv = J.inverse();
            const auto detJ = J.determinant();

            // Derivatives of the shape functions at the gauss node with respect to global coordinates x,y and z
            const Matrix<NumberOfNodes, 3, Eigen::RowMajor> dN_dx = (Jinv.transpose() * hexa.dL(
                gauss_point.first).transpose()).transpose();

            hexa_quadrature_nodes.push_back(GaussNode({
                                                          gauss_point.second,
                                                          detJ,
                                                          dN_dx,
                                                          Mat33::Identity()
                                                      }));
        }
    }

    Real v = 0.;
    for (std::size_t hexa_id = 0; hexa_id < topology->getNbHexahedra(); ++hexa_id) {
        for (std::size_t gauss_node_id = 0; gauss_node_id < p_quadrature_nodes[hexa_id].size(); ++gauss_node_id) {
            v += p_quadrature_nodes[hexa_id][gauss_node_id].weight*p_quadrature_nodes[hexa_id][gauss_node_id].jacobian_determinant;
        }
    }
    msg_info() << "Total volume is " << v;

    // Initialize the stiffness matrix of every hexahedrons
    p_stiffness_matrices.resize(topology->getNbHexahedra());
    p_initial_rotation.resize(topology->getNbHexahedra(), Mat33::Identity());
    p_current_rotation.resize(topology->getNbHexahedra(), Mat33::Identity());

    // Initialize the initial frame of each hexahedron
    if (d_corotated.getValue()) {
        for (std::size_t hexa_id = 0; hexa_id < topology->getNbHexahedra(); ++hexa_id) {
            Hexahedron hexa = hexahedron(hexa_id, X);
            p_initial_rotation[hexa_id] = hexa.frame({0, 0, 0});
        }
    }

    // Compute the initial tangent stiffness matrix
    compute_K();
}

void HexahedronElasticForce::addForce(
        const MechanicalParams * mparams,
        Data<VecDeriv> & d_f,
        const Data<VecCoord> & d_x,
        const Data<VecDeriv> & d_v)
{
    SOFA_UNUSED(mparams);
    SOFA_UNUSED(d_v);

    auto topology = d_topology_container.get();
    MechanicalState<DataTypes> * state = this->mstate.get();

    if (!topology or !state)
        return;

    if (p_stiffness_matrices.size() != topology->getNbHexahedra())
        return;

    sofa::helper::ReadAccessor<Data<VecCoord>> x = d_x;
    sofa::helper::ReadAccessor<Data<VecCoord>> x0 =  state->readRestPositions();
    sofa::helper::WriteAccessor<Data<VecDeriv>> f = d_f;

    const std::vector<Rotation> & initial_rotation = p_initial_rotation;
    std::vector<Rotation> & current_rotation = p_current_rotation;

    bool corotated = d_corotated.getValue();

    sofa::helper::AdvancedTimer::stepBegin("HexahedronElasticForce::addForce");
    const auto number_of_elements = topology->getNbHexahedra();
    for (std::size_t hexa_id = 0; hexa_id < number_of_elements; ++hexa_id) {
        Hexahedron hexa = hexahedron(hexa_id, x);

        const Rotation & R0 = initial_rotation[hexa_id];
        const Rotation R0t = R0.transpose();

        Rotation & R = current_rotation[hexa_id];

        // Extract the hexahedron's frame
        if (corotated)
            R = hexa.frame({0, 0, 0});

        const Rotation Rt = R.transpose();

        // Gather the displacement vector
        Vec24 U;
        Eigen::Index i = 0;
        for (const auto &node_id : topology->getHexahedron(static_cast<Topology::HexaID>(hexa_id))) {
            const Vec3 r0 {x0[node_id][0], x0[node_id][1],  x0[node_id][2]};
            const Vec3 r  {x [node_id][0],  x [node_id][1], x [node_id][2]};

            const Vec3 u = Rt*r - R0t*r0;

            U[i++] = u[0];
            U[i++] = u[1];
            U[i++] = u[2];
        }

        // Compute the force vector
        const auto &K = p_stiffness_matrices[hexa_id];
        Vec24 F = K.template selfadjointView<Eigen::Upper>() * U;

        // Write the forces into the output vector
        i = 0;
        for (const auto &node_id : topology->getHexahedron(static_cast<Topology::HexaID>(hexa_id))) {
            Vec3 force {F[i*3+0], F[i*3+1], F[i*3+2]};
            force = R*force;

            f[node_id][0] -= force[0];
            f[node_id][1] -= force[1];
            f[node_id][2] -= force[2];
            ++i;
        }
    }
    sofa::helper::AdvancedTimer::stepEnd("HexahedronElasticForce::addForce");
}

void HexahedronElasticForce::addDForce(
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

    auto kFactor = static_cast<Real>(mparams->kFactorIncludingRayleighDamping(this->rayleighStiffness.getValue()));
    sofa::helper::ReadAccessor<Data<VecDeriv>> dx = d_dx;
    sofa::helper::WriteAccessor<Data<VecDeriv>> df = d_df;
    std::vector<Rotation> & current_rotation = p_current_rotation;

    sofa::helper::AdvancedTimer::stepBegin("HexahedronElasticForce::addDForce");
    const auto number_of_elements = topology->getNbHexahedra();
    for (std::size_t hexa_id = 0; hexa_id < number_of_elements; ++hexa_id) {

        const Rotation & R  = current_rotation[hexa_id];
        const Rotation   Rt = R.transpose();

        // Gather the displacement vector
        Vec24 U;
        Eigen::Index i = 0;
        for (const auto & node_id : topology->getHexahedron(static_cast<Topology::HexaID>(hexa_id))) {
            const Vec3 v = {dx[node_id][0], dx[node_id][1], dx[node_id][2]};
            const Vec3 u = Rt*v;

            U[i++] = u[0];
            U[i++] = u[1];
            U[i++] = u[2];
        }

        // Compute the force vector
        const auto & K = p_stiffness_matrices[hexa_id];
        Vec24 F = K.template selfadjointView<Eigen::Upper>()*U*kFactor;

        // Write the forces into the output vector
        i = 0;
        for (const auto & node_id : topology->getHexahedron(static_cast<Topology::HexaID>(hexa_id))) {
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

void HexahedronElasticForce::addKToMatrix(sofa::defaulttype::BaseMatrix * matrix, SReal kFact, unsigned int & offset)
{
    auto * topology = d_topology_container.get();

    if (!topology)
        return;

    std::vector<Rotation> & current_rotation = p_current_rotation;

    sofa::helper::AdvancedTimer::stepBegin("HexahedronElasticForce::addKToMatrix");

    const auto number_of_elements = topology->getNbHexahedra();
    for (std::size_t hexa_id = 0; hexa_id < number_of_elements; ++hexa_id) {
        const auto & node_indices = topology->getHexahedron(static_cast<Topology::HexaID>(hexa_id));
        sofa::defaulttype::Mat3x3 R;
        for (size_t m = 0; m < 3; ++m) {
            for (size_t n = 0; n < 3; ++n) {
                R(m, n) = current_rotation[hexa_id](m, n);
            }
        }
        const auto   Rt = R.transposed();

        // Since the matrix K is block symmetric, we only kept the DxD blocks on the upper-triangle the matrix.
        // Here we need to accumulate the full matrix into Sofa's BaseMatrix.
        const auto & K = p_stiffness_matrices[hexa_id];

        // Blocks on the diagonal
        for (size_t i = 0; i < NumberOfNodes; ++i) {
            const auto x = (i*3);
            sofa::defaulttype::Mat<3, 3, Real> Kii;
            for (size_t m = 0; m < 3; ++m) {
                for (size_t n = 0; n < 3; ++n) {
                    Kii(m, n) = K(x+m, x+n);
                }
            }

            Kii = -1. * R*Kii*Rt *kFact;

            matrix->add(offset+node_indices[i]*3, offset+node_indices[i]*3, Kii);
        }

        // Blocks on the upper triangle
        for (size_t i = 0; i < NumberOfNodes; ++i) {
            for (size_t j = i+1; j < NumberOfNodes; ++j) {
                const auto x = (i*3);
                const auto y = (j*3);

                sofa::defaulttype::Mat<3, 3, Real> Kij;
                for (size_t m = 0; m < 3; ++m) {
                    for (size_t n = 0; n < 3; ++n) {
                        Kij(m, n) = K(x+m, y+n);
                    }
                }

                Kij = -1. * R*Kij*Rt *kFact;

                matrix->add(offset+node_indices[i]*3, offset+node_indices[j]*3, Kij);
                matrix->add(offset+node_indices[j]*3, offset+node_indices[i]*3, Kij.transposed());
            }
        }
    }
    sofa::helper::AdvancedTimer::stepEnd("HexahedronElasticForce::addKToMatrix");
}

void HexahedronElasticForce::compute_K()
{
    auto topology = d_topology_container.get();

    if (!topology)
        return;

    if (p_stiffness_matrices.size() != topology->getNbHexahedra())
        return;

    const Real youngModulus = d_youngModulus.getValue();
    const Real poissonRatio = d_poissonRatio.getValue();

    const Real l = youngModulus * poissonRatio / ((1 + poissonRatio) * (1 - 2 * poissonRatio));
    const Real mu = youngModulus / (2 * (1 + poissonRatio));

    Eigen::Matrix<Real, 6, 6> C;
    C <<
      l + 2*mu,    l,          l,       0,  0,  0,
        l,       l + 2*mu,     l,       0,  0,  0,
        l,         l,        l + 2*mu,  0,  0,  0,
        0,         0,          0,      mu,  0,  0,
        0,         0,          0,       0, mu,  0,
        0,         0,          0,       0,  0, mu;
    sofa::helper::AdvancedTimer::stepBegin("HexahedronElasticForce::compute_k");

    const auto number_of_elements = topology->getNbHexahedra();
    for (std::size_t hexa_id = 0; hexa_id < number_of_elements; ++hexa_id) {
        auto & K = p_stiffness_matrices[hexa_id];
        K.fill(0.);

        for (GaussNode &gauss_node : p_quadrature_nodes[hexa_id]) {
            // Jacobian of the gauss node's transformation mapping from the elementary space to the world space
            const Real detJ = gauss_node.jacobian_determinant;

            // Derivatives of the shape functions at the gauss node with respect to global coordinates x,y and z
            const auto dN_dx = gauss_node.dN_dx;

            // Gauss quadrature node weight
            const auto w = gauss_node.weight;

            // Computation of the element tangent-stiffness matrix
            for (std::size_t i = 0; i < NumberOfNodes; ++i) {
                // Derivatives of the ith shape function at the gauss node with respect to global coordinates x,y and z
                const Vec3 dxi = dN_dx.row(i).transpose();

                Matrix<6,3> Bi;
                Bi <<
                    dxi[0],    0  ,    0  ,
                       0  , dxi[1],    0  ,
                       0  ,    0  , dxi[2],
                    dxi[1], dxi[0],    0  ,
                       0  , dxi[2], dxi[1],
                    dxi[2],    0  , dxi[0];
                for (std::size_t j = i; j < NumberOfNodes; ++j) {
                    // Derivatives of the jth shape function at the gauss node with respect to global coordinates x,y and z
                    const Vec3 dxj = dN_dx.row(j).transpose();
                    Matrix<6,3> Bj;
                    Bj <<
                        dxj[0],    0  ,    0  ,
                           0  , dxj[1],    0  ,
                           0  ,    0  , dxj[2],
                        dxj[1], dxj[0],    0  ,
                           0  , dxj[2], dxj[1],
                        dxj[2],    0  , dxj[0];

                    K.template block<3, 3>(i*3, j*3).noalias() += (Bi.transpose()*C*Bj) * detJ * w;
                }
            }
        }
    }
    K_is_up_to_date = false;
    eigenvalues_are_up_to_date = false;
    sofa::helper::AdvancedTimer::stepEnd("HexahedronElasticForce::compute_k");
}

const Eigen::SparseMatrix<HexahedronElasticForce::Real> & HexahedronElasticForce::K() {
    if (not K_is_up_to_date) {
        const sofa::helper::ReadAccessor<Data<VecCoord>> X = this->mstate->readRestPositions();
        const auto nDofs = X.size() * 3;
        p_K.resize(nDofs, nDofs);
        p_K.setZero();

        ///< Triplets are used to store matrix entries before the call to 'compress'.
        /// Duplicates entries are summed up.
        std::vector<Eigen::Triplet<Real >> triplets;
        triplets.reserve(nDofs*24*2);

        auto *topology = d_topology_container.get();

        if (topology) {

            const std::vector<Rotation> &current_rotation = p_current_rotation;

            const auto number_of_elements = topology->getNbHexahedra();
            for (std::size_t hexa_id = 0; hexa_id < number_of_elements; ++hexa_id) {
                const auto &node_indices = topology->getHexahedron(hexa_id);
                const Rotation &R = current_rotation[hexa_id];
                const Rotation Rt = R.transpose();

                // Since the matrix K is block symmetric, we only kept the DxD blocks on the upper-triangle the matrix.
                // Here we need to accumulate the full matrix into Sofa's BaseMatrix.
                const auto &Ke = p_stiffness_matrices[hexa_id];

                // Blocks on the diagonal
                for (size_t i = 0; i < NumberOfNodes; ++i) {
                    const auto x = (node_indices[i]*3);
                    const Mat33 Kii = -1 * R * Ke.block<3,3>(i,i) * Rt;
                    for (size_t m = 0; m < 3; ++m) {
                        for (size_t n = 0; n < 3; ++n) {
                            triplets.emplace_back(x+m, x+n, Kii(m,n));
                        }
                    }
                }

                // Blocks on the upper triangle
                for (size_t i = 0; i < NumberOfNodes; ++i) {
                    for (size_t j = i+1; j < NumberOfNodes; ++j) {
                        const auto x = (node_indices[i]*3);
                        const auto y = (node_indices[j]*3);

                        const Mat33 Kij = -1 * R * Ke.block<3,3>(i,j) * Rt;

                        for (size_t m = 0; m < 3; ++m) {
                            for (size_t n = 0; n < 3; ++n) {
                                triplets.emplace_back(x+m, y+n, Kij(m,n));
                            }
                        }
                    }
                }
            }
        }
        p_K.setFromTriplets(triplets.begin(), triplets.end());
        K_is_up_to_date = true;
    }

    return p_K;
}

const Eigen::Matrix<HexahedronElasticForce::Real, Eigen::Dynamic, 1> & HexahedronElasticForce::eigenvalues()
{
    if (not eigenvalues_are_up_to_date) {
#ifdef EIGEN_USE_LAPACKE
        Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> k (K());
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>> eigensolver(k, Eigen::EigenvaluesOnly);
#else
        Eigen::SelfAdjointEigenSolver<Eigen::SparseMatrix<Real>> eigensolver(K(), Eigen::EigenvaluesOnly);
#endif
        if (eigensolver.info() != Eigen::Success) {
            msg_error() << "Unable to find the eigen values of K.";
        }

        p_eigenvalues = eigensolver.eigenvalues();
        eigenvalues_are_up_to_date = true;
    }

    return p_eigenvalues;
}

HexahedronElasticForce::Real HexahedronElasticForce::cond()
{
    const auto & values = eigenvalues();
    const auto min = values.minCoeff();
    const auto max = values.maxCoeff();

    return min/max;
}

void HexahedronElasticForce::computeBBox(const sofa::core::ExecParams*, bool onlyVisible)
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

    this->f_bbox.setValue(sofa::defaulttype::TBoundingBox<Real>(minBBox,maxBBox));
}

void HexahedronElasticForce::draw(const sofa::core::visual::VisualParams* vparams)
{
    auto * topology = d_topology_container.get();
    if (!topology)
        return;

    if (!vparams->displayFlags().getShowForceFields())
        return;

    vparams->drawTool()->saveLastState();

    if (vparams->displayFlags().getShowWireFrame())
        vparams->drawTool()->setPolygonMode(0,true);

    vparams->drawTool()->disableLighting();

    const VecCoord& x = this->mstate->read(sofa::core::ConstVecCoordId::position())->getValue();

    std::vector< sofa::defaulttype::Vector3 > points[6];
    const auto number_of_elements = topology->getNbHexahedra();
    for (std::size_t hexa_id = 0; hexa_id < number_of_elements; ++hexa_id) {
        const auto & node_indices = topology->getHexahedron(hexa_id);

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
    }

    vparams->drawTool()->drawTriangles(points[0], sofa::defaulttype::Vec<4,float>(0.7f,0.7f,0.1f,1.0f));
    vparams->drawTool()->drawTriangles(points[1], sofa::defaulttype::Vec<4,float>(0.7f,0.0f,0.0f,1.0f));
    vparams->drawTool()->drawTriangles(points[2], sofa::defaulttype::Vec<4,float>(0.0f,0.7f,0.0f,1.0f));
    vparams->drawTool()->drawTriangles(points[3], sofa::defaulttype::Vec<4,float>(0.0f,0.0f,0.7f,1.0f));
    vparams->drawTool()->drawTriangles(points[4], sofa::defaulttype::Vec<4,float>(0.1f,0.7f,0.7f,1.0f));
    vparams->drawTool()->drawTriangles(points[5], sofa::defaulttype::Vec<4,float>(0.7f,0.1f,0.7f,1.0f));


    std::vector< sofa::defaulttype::Vector3 > ignored_points[6];
    for (std::size_t hexa_id = 0; hexa_id < number_of_elements; ++hexa_id) {
        const auto & node_indices = topology->getHexahedron(hexa_id);

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



        ignored_points[0].push_back(pa);
        ignored_points[0].push_back(pb);
        ignored_points[0].push_back(pc);
        ignored_points[0].push_back(pa);
        ignored_points[0].push_back(pc);
        ignored_points[0].push_back(pd);

        ignored_points[1].push_back(pe);
        ignored_points[1].push_back(pf);
        ignored_points[1].push_back(pg);
        ignored_points[1].push_back(pe);
        ignored_points[1].push_back(pg);
        ignored_points[1].push_back(ph);

        ignored_points[2].push_back(pc);
        ignored_points[2].push_back(pd);
        ignored_points[2].push_back(ph);
        ignored_points[2].push_back(pc);
        ignored_points[2].push_back(ph);
        ignored_points[2].push_back(pg);

        ignored_points[3].push_back(pa);
        ignored_points[3].push_back(pb);
        ignored_points[3].push_back(pf);
        ignored_points[3].push_back(pa);
        ignored_points[3].push_back(pf);
        ignored_points[3].push_back(pe);

        ignored_points[4].push_back(pa);
        ignored_points[4].push_back(pd);
        ignored_points[4].push_back(ph);
        ignored_points[4].push_back(pa);
        ignored_points[4].push_back(ph);
        ignored_points[4].push_back(pe);

        ignored_points[5].push_back(pb);
        ignored_points[5].push_back(pc);
        ignored_points[5].push_back(pg);
        ignored_points[5].push_back(pb);
        ignored_points[5].push_back(pg);
        ignored_points[5].push_back(pf);
    }

    vparams->drawTool()->drawTriangles(ignored_points[0], sofa::defaulttype::Vec<4,float>(0.49f,0.49f,0.49f,0.3f));
    vparams->drawTool()->drawTriangles(ignored_points[1], sofa::defaulttype::Vec<4,float>(0.49f,0.49f,0.49f,0.3f));
    vparams->drawTool()->drawTriangles(ignored_points[2], sofa::defaulttype::Vec<4,float>(0.49f,0.49f,0.49f,0.3f));
    vparams->drawTool()->drawTriangles(ignored_points[3], sofa::defaulttype::Vec<4,float>(0.49f,0.49f,0.49f,0.3f));
    vparams->drawTool()->drawTriangles(ignored_points[4], sofa::defaulttype::Vec<4,float>(0.49f,0.49f,0.49f,0.3f));
    vparams->drawTool()->drawTriangles(ignored_points[5], sofa::defaulttype::Vec<4,float>(0.49f,0.49f,0.49f,0.3f));


    if (vparams->displayFlags().getShowWireFrame())
        vparams->drawTool()->setPolygonMode(0,false);

    vparams->drawTool()->restoreLastState();
}


static int HexahedronElasticForceClass = RegisterObject("Caribou Hexahedron FEM Forcefield")
    .add< HexahedronElasticForce >(true)
;

} // namespace SofaCaribou::forcefield
