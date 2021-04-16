#include <numeric>
#include <queue>
#include <array>
#ifdef CARIBOU_WITH_OPENMP
#include <omp.h>
#endif

#include "FictitiousGridElasticForce.h"

DISABLE_ALL_WARNINGS_BEGIN
#include <sofa/core/visual/VisualParams.h>
#include <sofa/core/ObjectFactory.h>
#include <sofa/core/MechanicalParams.h>
#include <sofa/simulation/Node.h>
#include <sofa/helper/AdvancedTimer.h>
DISABLE_ALL_WARNINGS_END

#include <Eigen/Sparse>


#include <Caribou/Mechanics/Elasticity/Strain.h>



namespace SofaCaribou::forcefield {

using namespace sofa::core::topology;
using namespace caribou::geometry;
using namespace caribou::mechanics;


FictitiousGridElasticForce::FictitiousGridElasticForce()
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
    , d_integration_method(initData(&d_integration_method,
                                    "integration_method",
                                    R"(
                Integration method used to integrate the stiffness matrix.

                Methods are:
                  Regular:          Regular 8 points gauss integration.
                  SubdividedVolume: Hexas are recursively subdivided into cubic subcells and these subcells are used to
                  compute the inside volume of the regular hexa's gauss points.
                  ** Requires a sparse grid topology **
                  SubdividedGauss:  Hexas are recursively subdivided into cubic subcells and these subcells are used to
                  add new gauss points. Gauss points outside of the boundary are ignored.
                  ** Requires a sparse grid topology **
                )",
                                    true /*displayed_in_GUI*/, false /*read_only_in_GUI*/))
    , d_grid_container(initLink(
        "fictitious_grid", "Fictitious grid that contains the elements on which this force will be computed."))
{
    d_integration_method.setValue(sofa::helper::OptionsGroup(std::vector<std::string> {
        "Regular", "SubdividedVolume", "SubdividedGauss"
    }));

    sofa::helper::WriteAccessor<Data< sofa::helper::OptionsGroup >> integration_method = d_integration_method;
    integration_method->setSelectedItem((unsigned int) 0);
}

void FictitiousGridElasticForce::init()
{
    Inherit::init();
    if (not d_grid_container.get()) {
        auto containers = this->getContext()->template getObjects<FictitiousGrid>(BaseContext::Local);
        auto node = dynamic_cast<const sofa::simulation::Node *> (this->getContext());
        if (containers.empty()) {
            msg_error() << "No fictitious grid were found in the context node '" << node->getPathName() << "'.";
        } else if (containers.size() > 1) {
            msg_error() <<
                        "Multiple fictitious grids were found in the node '" << node->getPathName() << "'." <<
                        " Please specify which one contains the elements on which this force field will be computed " <<
                        "by explicitly setting the container's path in the 'fictitious_grid' parameter.";
        } else {
            d_grid_container.set(containers[0]);
            msg_info() << "Automatically found the fictitious grid at '" << d_grid_container.get()->getPathName() << "'.";
        }
    }

    if (!this->mstate.get())
        msg_error() << "No mechanical state set. Please add a mechanical object in the current node or specify a path "
                       "to one in the 'mstate' parameter.";


    reinit();
}

void FictitiousGridElasticForce::reinit()
{
    FictitiousGrid * grid = d_grid_container.get();
    MechanicalState<DataTypes> * state = this->mstate.get();

    if (!grid or !state)
        return;

    if (grid->number_of_cells() == 0) {
        msg_warning() << "The fictitious grid ('" << grid->getPathName() << "') does not contain any hexahedron.";
        return;
    }

    if (integration_method() == IntegrationMethod::SubdividedVolume or
        integration_method() == IntegrationMethod::SubdividedGauss) {

        if (grid->number_of_subdivisions() == 0) {
            msg_error() << "Integration method '" << integration_method_as_string()
                        << "' requires the fictitious grid to have cell subdivisions.";
            set_integration_method(IntegrationMethod::Regular);
        }
    }

    const sofa::helper::ReadAccessor<Data<VecCoord>> X = state->readRestPositions();

    // Make sure every node of the hexahedrons have its coordinates inside the mechanical state vector
    for (std::size_t hexa_id = 0; hexa_id < grid->number_of_cells(); ++hexa_id) {
        const auto & node_indices = grid->get_node_indices_of(hexa_id);
        for (sofa::Index j = 0; j < 8; ++j) {
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
    p_quadrature_nodes.resize(grid->number_of_cells());
    const auto int_method = integration_method();
    for (std::size_t hexa_id = 0; hexa_id < grid->number_of_cells(); ++hexa_id) {
        const auto e = grid->get_cell_element(hexa_id);

        if (int_method == IntegrationMethod::Regular) {
            p_quadrature_nodes[hexa_id].resize(caribou::geometry::traits<Hexahedron>::NumberOfGaussNodesAtCompileTime);
            for (std::size_t gauss_node_id = 0; gauss_node_id < caribou::geometry::traits<Hexahedron>::NumberOfGaussNodesAtCompileTime; ++gauss_node_id) {
                const auto & g = e.gauss_node(gauss_node_id);

                const auto J = e.jacobian(g.position);
                const Mat33 Jinv = J.inverse();
                const auto detJ = J.determinant();

                p_quadrature_nodes[hexa_id][gauss_node_id].weight = detJ * g.weight;
                p_quadrature_nodes[hexa_id][gauss_node_id].dN_dx = (Jinv.transpose() * e.dL(g.position).transpose()).transpose();
            }
        } else {
            const auto level = (int_method == IntegrationMethod::SubdividedVolume) ? 0 : grid->number_of_subdivisions();
            const auto gauss_nodes = grid->get_gauss_nodes_of_cell(hexa_id, level);
            p_quadrature_nodes[hexa_id].resize(gauss_nodes.size());
            for (std::size_t gauss_node_id = 0; gauss_node_id < gauss_nodes.size(); ++gauss_node_id) {
                const auto & gauss_node = gauss_nodes[gauss_node_id].first;
                const auto & gauss_weight = gauss_nodes[gauss_node_id].second;
                const auto J = e.jacobian(gauss_node);
                const Mat33 Jinv = J.inverse();

                p_quadrature_nodes[hexa_id][gauss_node_id].weight = gauss_weight;
                p_quadrature_nodes[hexa_id][gauss_node_id].dN_dx = (Jinv.transpose() * e.dL(gauss_node).transpose()).transpose();
            }
        }
    }

    Real v = 0.;
    UNSIGNED_INTEGER_TYPE negative_jacobians = 0;
    for (std::size_t hexa_id = 0; hexa_id < grid->number_of_cells(); ++hexa_id) {
        for (auto & gauss_node : p_quadrature_nodes[hexa_id]) {
            v += gauss_node.weight;
            
            if (gauss_node.weight < 0)
                negative_jacobians++;
        }
    }
    msg_info() << "Total volume is " << v;
    
    if (negative_jacobians > 0) {
        msg_warning() << negative_jacobians << " gauss points have a negative jacobian";
    }


    // Initialize the stiffness matrix of every hexahedrons
    p_stiffness_matrices.resize(grid->number_of_cells());
    p_initial_rotation.resize(grid->number_of_cells(), Mat33::Identity());
    p_current_rotation.resize(grid->number_of_cells(), Mat33::Identity());

    // Initialize the initial frame of each hexahedron
    if (d_corotated.getValue()) {
        if (not d_linear_strain.getValue()) {
            msg_warning() << "The corotated method won't be computed since nonlinear strain is used.";
        } else {
            for (std::size_t hexa_id = 0; hexa_id < grid->number_of_cells(); ++hexa_id) {
                Hexahedron hexa = hexahedron(hexa_id, X);
                p_initial_rotation[hexa_id] = hexa.frame({0, 0, 0});
            }
        }
    }

#ifdef CARIBOU_WITH_OPENMP
    const auto * env = std::getenv("OMP_NUM_THREADS");
    if (env) {
        msg_info() << "Using " << env << " threads for computations.";
    } else {
        msg_info() << "Using 1 threads for computations.";
    }
#endif

    // Compute the initial tangent stiffness matrix
    compute_K();
}

void FictitiousGridElasticForce::addForce(
    const MechanicalParams * mparams,
    Data<VecDeriv> & d_f,
    const Data<VecCoord> & d_x,
    const Data<VecDeriv> & d_v)
{
    SOFA_UNUSED(mparams);
    SOFA_UNUSED(d_v);

    Mat33 Id = Mat33::Identity();

    auto grid = d_grid_container.get();
    MechanicalState<DataTypes> * state = this->mstate.get();

    if (!grid or !state)
        return;

    if (p_stiffness_matrices.size() != grid->number_of_cells())
        return;

    sofa::helper::ReadAccessor<Data<VecCoord>> x = d_x;
    sofa::helper::ReadAccessor<Data<VecCoord>> x0 =  state->readRestPositions();
    sofa::helper::WriteAccessor<Data<VecDeriv>> f = d_f;

    const std::vector<Mat33> & initial_rotation = p_initial_rotation;
    std::vector<Mat33> & current_rotation = p_current_rotation;

    bool corotated = d_corotated.getValue();
    bool linear = d_linear_strain.getValue();
    recompute_compute_tangent_stiffness = (not linear);
    if (linear) {
        // Small (linear) strain
        sofa::helper::AdvancedTimer::stepBegin("FictitiousGridElasticForce::addForce");
#pragma omp parallel for
        for (int hexa_id = 0; hexa_id < static_cast<int>(grid->number_of_cells()); ++hexa_id) {
            Hexahedron hexa = hexahedron(hexa_id, x);

            const Mat33 & R0 = initial_rotation[hexa_id];
            const Mat33   R0t = R0.transpose();

            Mat33 & R = current_rotation[hexa_id];


            // Extract the hexahedron's frame
            if (corotated)
                R = hexa.frame({0, 0, 0});

            const Mat33 Rt = R.transpose();

            // Gather the displacement vector
            Vec24 U;
            size_t i = 0;
            for (const auto &node_id : grid->get_node_indices_of(hexa_id)) {
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
            for (const auto &node_id : grid->get_node_indices_of(hexa_id)) {
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
        sofa::helper::AdvancedTimer::stepEnd("FictitiousGridElasticForce::addForce");
    } else {
        // Nonlinear Green-Lagrange strain
        sofa::helper::AdvancedTimer::stepBegin("FictitiousGridElasticForce::addForce");
        const Real youngModulus = d_youngModulus.getValue();
        const Real poissonRatio = d_poissonRatio.getValue();

        const Real l = youngModulus * poissonRatio / ((1 + poissonRatio) * (1 - 2 * poissonRatio));
        const Real m = youngModulus / (2 * (1 + poissonRatio));

#pragma omp parallel for
        for (int hexa_id = 0; hexa_id < static_cast<int>(grid->number_of_cells()); ++hexa_id) {
            const auto &hexa = grid->get_node_indices_of(hexa_id);

            Matrix<8, 3, Eigen::RowMajor> U;

            for (sofa::Index i = 0; i < 8; ++i) {
                const auto u = x[hexa[i]] - x0[hexa[i]];
                U(i, 0) = u[0];
                U(i, 1) = u[1];
                U(i, 2) = u[2];
            }

            Matrix<8, 3, Eigen::RowMajor> forces;
            forces.fill(0);
            for (GaussNode &gauss_node : p_quadrature_nodes[hexa_id]) {
                // Derivatives of the shape functions at the gauss node with respect to global coordinates x,y and z
                const auto dN_dx = gauss_node.dN_dx;

                // Gauss quadrature node weight
                const auto w = gauss_node.weight;

                // Deformation tensor at gauss node
                auto & F = gauss_node.F;
                F = elasticity::strain::F(dN_dx, U);

                // Strain tensor at gauss node
                const Mat33 C = F.transpose()*F;
                const Mat33 E = 1/2. * (C - Id);

                // Stress tensor at gauss node
                const Mat33 S = 2.*m*E + (l * E.trace() * Id);

                // Elastic forces w.r.t the gauss node applied on each nodes
                for (Eigen::Index i = 0; i < 8; ++i) {
                    const Vec3 dx = dN_dx.row(i).transpose();
                    const Vec3 f_ = w * (F*S) * dx;
                    forces(i, 0) += f_[0];
                    forces(i, 1) += f_[1];
                    forces(i, 2) += f_[2];
                }
            }

            for (Eigen::Index i = 0; i < 8; ++i) {
#pragma omp atomic
                f[hexa[static_cast<sofa::Index>(i)]][0] -= forces.row(i)[0];

#pragma omp atomic
                f[hexa[static_cast<sofa::Index>(i)]][1] -= forces.row(i)[1];

#pragma omp atomic
                f[hexa[static_cast<sofa::Index>(i)]][2] -= forces.row(i)[2];
            }
        }
        sofa::helper::AdvancedTimer::stepEnd("FictitiousGridElasticForce::addForce");
    }
}

void FictitiousGridElasticForce::addDForce(
    const MechanicalParams* mparams,
    Data<VecDeriv>& d_df,
    const Data<VecDeriv>& d_dx)
{
    auto * grid = d_grid_container.get();
    MechanicalState<DataTypes> * state = this->mstate.get();

    if (!grid or !state)
        return;

    if (p_stiffness_matrices.size() != grid->number_of_cells())
        return;

    if (recompute_compute_tangent_stiffness)
        compute_K();

    auto kFactor = (Real)mparams->kFactorIncludingRayleighDamping(this->rayleighStiffness.getValue());
    sofa::helper::ReadAccessor<Data<VecDeriv>> dx = d_dx;
    sofa::helper::WriteAccessor<Data<VecDeriv>> df = d_df;
    std::vector<Mat33> & current_rotation = p_current_rotation;

    sofa::helper::AdvancedTimer::stepBegin("FictitiousGridElasticForce::addDForce");
#pragma omp parallel for
    for (int hexa_id = 0; hexa_id < static_cast<int>(grid->number_of_cells()); ++hexa_id) {
        const Mat33 & R  = current_rotation[hexa_id];
        const Mat33   Rt = R.transpose();

        // Gather the displacement vector
        Vec24 U;
        size_t i = 0;
        for (const auto & node_id : grid->get_node_indices_of(hexa_id)) {
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
        for (const auto & node_id : grid->get_node_indices_of(hexa_id)) {
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
    sofa::helper::AdvancedTimer::stepEnd("FictitiousGridElasticForce::addDForce");
}

void FictitiousGridElasticForce::addKToMatrix(sofa::defaulttype::BaseMatrix * matrix, SReal kFact, unsigned int & offset)
{
    auto * grid = d_grid_container.get();

    if (!grid)
        return;

    if (recompute_compute_tangent_stiffness)
        compute_K();

    std::vector<Mat33> & current_rotation = p_current_rotation;

    sofa::helper::AdvancedTimer::stepBegin("FictitiousGridElasticForce::addKToMatrix");

    const auto number_of_elements = grid->number_of_cells();
#pragma omp parallel for
    for (int hexa_id = 0; hexa_id < static_cast<int>(number_of_elements); ++hexa_id) {
        const auto & node_indices = grid->get_node_indices_of(hexa_id);
        sofa::defaulttype::Mat3x3 R;
        for (unsigned int m = 0; m < 3; ++m) {
            for (unsigned int n = 0; n < 3; ++n) {
                R(m, n) = current_rotation[hexa_id](m, n);
            }
        }
        const auto   Rt = R.transposed();

        // Since the matrix K is block symmetric, we only kept the DxD blocks on the upper-triangle the matrix.
        // Here we need to accumulate the full matrix into Sofa's BaseMatrix.
        const auto & K = p_stiffness_matrices[hexa_id];

        // Blocks on the diagonal
        for (sofa::Index i = 0; i < NumberOfNodes; ++i) {
            const auto x = static_cast<unsigned int>(i*3);
            sofa::defaulttype::Mat<3, 3, Real> Kii;
            for (unsigned int m = 0; m < 3; ++m) {
                for (unsigned int n = 0; n < 3; ++n) {
                    Kii(m, n) = K(x+m, x+n);
                }
            }

            Kii = -1. * R*Kii*Rt *kFact;
#pragma omp critical
            matrix->add(offset+node_indices[i]*3, offset+node_indices[i]*3, Kii);
        }

        // Blocks on the upper triangle
        for (sofa::Index i = 0; i < NumberOfNodes; ++i) {
            for (sofa::Index j = i+1; j < NumberOfNodes; ++j) {
                const auto x = static_cast<unsigned int>(i*3);
                const auto y = static_cast<unsigned int>(j*3);

                sofa::defaulttype::Mat<3, 3, Real> Kij;
                for (unsigned int m = 0; m < 3; ++m) {
                    for (unsigned int n = 0; n < 3; ++n) {
                        Kij(m, n) = K(x+m, y+n);
                    }
                }

                Kij = -1. * R*Kij*Rt *kFact;
#pragma omp critical
                matrix->add(offset+node_indices[i]*3, offset+node_indices[j]*3, Kij);
#pragma omp critical
                matrix->add(offset+node_indices[j]*3, offset+node_indices[i]*3, Kij.transposed());
            }
        }
    }
    sofa::helper::AdvancedTimer::stepEnd("FictitiousGridElasticForce::addKToMatrix");
}

void FictitiousGridElasticForce::compute_K()
{
    auto grid = d_grid_container.get();

    if (!grid)
        return;

    if (p_stiffness_matrices.size() != grid->number_of_cells())
        return;

    Mat33 Id = Mat33::Identity();
    const Real youngModulus = d_youngModulus.getValue();
    const Real poissonRatio = d_poissonRatio.getValue();

    const Real l = youngModulus * poissonRatio / ((1 + poissonRatio) * (1 - 2 * poissonRatio));
    const Real m = youngModulus / (2 * (1 + poissonRatio));

    Eigen::Matrix<Real, 6, 6> C;
    C <<
      l + 2*m,     l,          l,       0,  0,  0,
        l,       l + 2*m,      l,       0,  0,  0,
        l,         l,        l + 2*m,   0,  0,  0,
        0,         0,          0,       m,  0,  0,
        0,         0,          0,       0, m,  0,
        0,         0,          0,       0,  0, m;

    sofa::helper::AdvancedTimer::stepBegin("FictitiousGridElasticForce::compute_k");

#pragma omp parallel for
    for (int hexa_id = 0; hexa_id < static_cast<int>(grid->number_of_cells()); ++hexa_id) {
        auto & K = p_stiffness_matrices[hexa_id];
        K.fill(0.);

        for (GaussNode &gauss_node : p_quadrature_nodes[hexa_id]) {
            // Derivatives of the shape functions at the gauss node with respect to global coordinates x,y and z
            const auto dN_dx = gauss_node.dN_dx;

            // Gauss quadrature node weight
            const auto w = gauss_node.weight;

            // Deformation tensor at gauss node
            Mat33 & F = gauss_node.F;

            // Strain tensor at gauss node
            const Mat33 E = 1/2. * (F.transpose()*F - Id);

            // Stress tensor at gauss node
            const Mat33 S = 2.*m*E + (l * E.trace() * Id);

            // Computation of the tangent-stiffness matrix
            for (std::size_t i = 0; i < 8; ++i) {
                // Derivatives of the ith shape function at the gauss node with respect to global coordinates x,y and z
                const Vec3 dxi = dN_dx.row(i).transpose();

                Matrix<6,3> Bi;
                Bi <<
                    F(0,0)*dxi[0],                 F(1,0)*dxi[0],                 F(2,0)*dxi[0],
                    F(0,1)*dxi[1],                 F(1,1)*dxi[1],                 F(2,1)*dxi[1],
                    F(0,2)*dxi[2],                 F(1,2)*dxi[2],                 F(2,2)*dxi[2],
                    F(0,0)*dxi[1] + F(0,1)*dxi[0], F(1,0)*dxi[1] + F(1,1)*dxi[0], F(2,0)*dxi[1] + F(2,1)*dxi[0],
                    F(0,1)*dxi[2] + F(0,2)*dxi[1], F(1,1)*dxi[2] + F(1,2)*dxi[1], F(2,1)*dxi[2] + F(2,2)*dxi[1],
                    F(0,0)*dxi[2] + F(0,2)*dxi[0], F(0,1)*dxi[2] + F(1,2)*dxi[0], F(2,0)*dxi[2] + F(2,2)*dxi[0];

                for (std::size_t j = i; j < 8; ++j) {
                    // Derivatives of the jth shape function at the gauss node with respect to global coordinates x,y and z
                    const Vec3 dxj = dN_dx.row(j).transpose();

                    Matrix<6,3> Bj;
                    Bj <<
                        F(0,0)*dxj[0],                 F(1,0)*dxj[0],                 F(2,0)*dxj[0],
                        F(0,1)*dxj[1],                 F(1,1)*dxj[1],                 F(2,1)*dxj[1],
                        F(0,2)*dxj[2],                 F(1,2)*dxj[2],                 F(2,2)*dxj[2],
                        F(0,0)*dxj[1] + F(0,1)*dxj[0], F(1,0)*dxj[1] + F(1,1)*dxj[0], F(2,0)*dxj[1] + F(2,1)*dxj[0],
                        F(0,1)*dxj[2] + F(0,2)*dxj[1], F(1,1)*dxj[2] + F(1,2)*dxj[1], F(2,1)*dxj[2] + F(2,2)*dxj[1],
                        F(0,0)*dxj[2] + F(0,2)*dxj[0], F(0,1)*dxj[2] + F(1,2)*dxj[0], F(2,0)*dxj[2] + F(2,2)*dxj[0];

                    K.template block<3, 3>(i*3, j*3).noalias() += (dxi.dot(S*dxj)*Id + Bi.transpose()*C*Bj)  * w;
                }
            }
        }
    }
    recompute_compute_tangent_stiffness = false;
    K_is_up_to_date = false;
    eigenvalues_are_up_to_date = false;
    sofa::helper::AdvancedTimer::stepEnd("FictitiousGridElasticForce::compute_k");
}

const Eigen::SparseMatrix<FictitiousGridElasticForce::Real> & FictitiousGridElasticForce::K() {
    if (not K_is_up_to_date) {
        const sofa::helper::ReadAccessor<Data<VecCoord>> X = this->mstate->readRestPositions();
        const auto nDofs = X.size() * 3;
        p_K.resize(nDofs, nDofs);
        p_K.setZero();

        ///< Triplets are used to store matrix entries before the call to 'compress'.
        /// Duplicates entries are summed up.
        std::vector<Eigen::Triplet<Real >> triplets;
        triplets.reserve(nDofs*24*2);

        auto *grid = d_grid_container.get();

        if (grid) {

            const std::vector<Mat33> &current_rotation = p_current_rotation;

            const auto number_of_elements = grid->number_of_cells();
            for (std::size_t hexa_id = 0; hexa_id < number_of_elements; ++hexa_id) {
                const auto & node_indices = grid->get_node_indices_of(hexa_id);
                const Mat33 &R = current_rotation[hexa_id];
                const Mat33 Rt = R.transpose();

                // Since the matrix K is block symmetric, we only kept the DxD blocks on the upper-triangle the matrix.
                // Here we need to accumulate the full matrix.
                const auto &Ke = p_stiffness_matrices[hexa_id];


                // Blocks on the diagonal
                for (sofa::Index i = 0; i < NumberOfNodes; ++i) {
                    const auto x = static_cast<int>(node_indices[i]*3);
                    const Mat33 Kii = -1 * R * Ke.block<3,3>(i,i) * Rt;
                    for (int m = 0; m < 3; ++m) {
                        for (int n = 0; n < 3; ++n) {
                            triplets.emplace_back(x+m, x+n, Kii(m,n));
                        }
                    }
                }

                // Blocks on the upper triangle
                for (sofa::Index i = 0; i < NumberOfNodes; ++i) {
                    for (sofa::Index j = i+1; j < NumberOfNodes; ++j) {
                        const auto x = static_cast<int>(node_indices[i]*3);
                        const auto y = static_cast<int>(node_indices[j]*3);

                        const Mat33 Kij = -1 * R * Ke.block<3,3>(i,j) * Rt;
                        const Mat33 Kji = Kij.transpose();

                        for (int m = 0; m < 3; ++m) {
                            for (int n = 0; n < 3; ++n) {
                                triplets.emplace_back(x+m, y+n, Kij(m,n));
                                triplets.emplace_back(y+m, x+n, Kji(m,n));
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

const Eigen::Matrix<FictitiousGridElasticForce::Real, Eigen::Dynamic, 1> & FictitiousGridElasticForce::eigenvalues()
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

FictitiousGridElasticForce::Real FictitiousGridElasticForce::cond()
{
    const auto & values = eigenvalues();
    const auto min = values.minCoeff();
    const auto max = values.maxCoeff();

    return min/max;
}

void FictitiousGridElasticForce::computeBBox(const sofa::core::ExecParams*, bool onlyVisible)
{
    if( !onlyVisible ) return;

    sofa::helper::ReadAccessor<Data<VecCoord>> x = this->mstate->read(sofa::core::VecCoordId::position());

    const Real max_real = std::numeric_limits<Real>::max();
    const Real min_real = std::numeric_limits<Real>::lowest();
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

void FictitiousGridElasticForce::draw(const sofa::core::visual::VisualParams* vparams)
{
    auto * grid = d_grid_container.get();
    if (!grid)
        return;

    if (!vparams->displayFlags().getShowForceFields())
        return;

    vparams->drawTool()->saveLastState();

    if (vparams->displayFlags().getShowWireFrame())
        vparams->drawTool()->setPolygonMode(0,true);

    vparams->drawTool()->disableLighting();

    const VecCoord& x = this->mstate->read(sofa::core::ConstVecCoordId::position())->getValue();

    std::vector< sofa::defaulttype::Vector3 > points[6];
    for (std::size_t hexa_id = 0; hexa_id < grid->number_of_cells(); ++hexa_id) {
        const auto & node_indices = grid->get_node_indices_of(hexa_id);

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

    vparams->drawTool()->drawTriangles(points[0], sofa::helper::types::RGBAColor(0.7f,0.7f,0.1f,1.0f));
    vparams->drawTool()->drawTriangles(points[1], sofa::helper::types::RGBAColor(0.7f,0.0f,0.0f,1.0f));
    vparams->drawTool()->drawTriangles(points[2], sofa::helper::types::RGBAColor(0.0f,0.7f,0.0f,1.0f));
    vparams->drawTool()->drawTriangles(points[3], sofa::helper::types::RGBAColor(0.0f,0.0f,0.7f,1.0f));
    vparams->drawTool()->drawTriangles(points[4], sofa::helper::types::RGBAColor(0.1f,0.7f,0.7f,1.0f));
    vparams->drawTool()->drawTriangles(points[5], sofa::helper::types::RGBAColor(0.7f,0.1f,0.7f,1.0f));

    if (vparams->displayFlags().getShowWireFrame())
        vparams->drawTool()->setPolygonMode(0,false);

    vparams->drawTool()->restoreLastState();
}


static int FictitiousGridElasticForceClass = RegisterObject("Caribou Fictitious grid FEM Forcefield")
    .add< FictitiousGridElasticForce >(true)
;

} // namespace SofaCaribou::forcefield
