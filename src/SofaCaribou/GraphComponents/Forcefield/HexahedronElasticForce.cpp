#include <numeric>
#include <queue>
#include <array>
#ifdef CARIBOU_WITH_OPENMP
#include <omp.h>
#endif
#include <Eigen/Sparse>

#include <sofa/core/visual/VisualParams.h>
#include <sofa/core/ObjectFactory.h>
#include <sofa/simulation/Node.h>
#include <sofa/helper/AdvancedTimer.h>

#include <Caribou/Geometry/Hexahedron.h>
#include <Caribou/Geometry/RectangularHexahedron.h>
#include <Caribou/Mechanics/Elasticity/Strain.h>

#include "HexahedronElasticForce.h"


namespace SofaCaribou::GraphComponents::forcefield {

using namespace sofa::core::topology;
using namespace caribou::geometry;
using namespace caribou::mechanics;
using sofa::component::topology::SparseGridTopology;

/**
 * Split an hexahedron into 8 subcells. Subcells have their node position relative to the center of the outer hexahedron.
 */
template<typename Hexahedron>
std::array<RectangularHexahedron<typename Hexahedron::CanonicalElement>, 8>
split_in_local_hexahedrons(const Hexahedron & h) {
    using LocalHexahedron = RectangularHexahedron<typename Hexahedron::CanonicalElement>;
    const auto hx = (h.node(1) - h.node(0)).norm();
    const auto hy = (h.node(3) - h.node(0)).norm();
    const auto hz = (h.node(4) - h.node(0)).norm();
    return {{
        LocalHexahedron({-0.5, -0.5, -0.5}, {hx/2., hy/2., hz/2.}),
        LocalHexahedron({+0.5, -0.5, -0.5}, {hx/2., hy/2., hz/2.}),
        LocalHexahedron({+0.5, +0.5, -0.5}, {hx/2., hy/2., hz/2.}),
        LocalHexahedron({-0.5, +0.5, -0.5}, {hx/2., hy/2., hz/2.}),
        LocalHexahedron({-0.5, -0.5, +0.5}, {hx/2., hy/2., hz/2.}),
        LocalHexahedron({+0.5, -0.5, +0.5}, {hx/2., hy/2., hz/2.}),
        LocalHexahedron({+0.5, +0.5, +0.5}, {hx/2., hy/2., hz/2.}),
        LocalHexahedron({-0.5, +0.5, +0.5}, {hx/2., hy/2., hz/2.})
    }};
}

/**
 * Recursively get all the gauss positions that are inside the surface domain with their correct weight.
 * @tparam Hexahedron
 * @param grid
 * @param h
 * @return An array of 8 vectors of gauss points (one per subcell). A gauss point is represented by the tuple (xi, w)
 *         where xi is the local coordinates vector of the gauss point and w is its weight.
 */
template<typename Hexahedron>
std::array<std::vector<std::pair<typename Hexahedron::LocalCoordinates, typename Hexahedron::Real>>, 8>
recursively_get_subcells_gauss_points(SparseGridTopology & grid, const Hexahedron & hexa, const std::size_t & max_number_of_subdivisions) {
    using LocalHexahedron = RectangularHexahedron<typename Hexahedron::CanonicalElement>;
    using LocalCoordinates = typename Hexahedron::LocalCoordinates;
    using Real = typename Hexahedron::Real;
    using MapVector = Eigen::Map<const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1>, Eigen::ColMajor>;

    std::array<std::vector<std::pair<typename Hexahedron::LocalCoordinates, typename Hexahedron::Real>>, 8> gauss_points;

    // First, check if all the top level hexa gauss points are inside to avoid the subdivision
    {
        bool all_local_gauss_points_are_inside = true;
        for (std::size_t gauss_node_id = 0; gauss_node_id < Hexahedron::number_of_gauss_nodes; ++gauss_node_id) {
            const auto &gauss_node = MapVector(Hexahedron::gauss_nodes[gauss_node_id]);
            const auto gauss_position = hexa.T(gauss_node);
            Real fx, fy, fz; // unused
            const auto gauss_cube_id = grid.findCube({gauss_position[0], gauss_position[1], gauss_position[2]}, fx, fy,fz);
            if (gauss_cube_id < 0) {
                all_local_gauss_points_are_inside = false;
            }
        }

        if (all_local_gauss_points_are_inside) {
            // They are all inside the boundary, let's add them and return since the hexa is not cut by the boundary
            for (std::size_t gauss_node_id = 0; gauss_node_id < Hexahedron::number_of_gauss_nodes; ++gauss_node_id) {
                const auto &gauss_node = MapVector(Hexahedron::gauss_nodes[gauss_node_id]);
                const auto &gauss_weight = Hexahedron::gauss_weights[gauss_node_id];
                gauss_points[gauss_node_id].push_back({gauss_node, gauss_weight});
            }
            return gauss_points;
        }
    }

    // At this point, we got an hexa that is either outside, or on the boundary of the surface
    const auto subcells = split_in_local_hexahedrons(LocalHexahedron({0,0,0}));

    for (std::size_t i = 0; i < subcells.size(); ++i) {

        // FIFO which will store final sub-cells. A sub-cell is final when it is fully inside
        // the surface domain, or when the maximum number of subdivisions has been reach. If it is not final,
        // then the sub-cell is removed from the top of the queue and its 8 sub-cells are added to the end
        // of the queue.
        struct SubCell {
            LocalHexahedron hexahedron;
            UNSIGNED_INTEGER_TYPE subdivision_level;
        };

        std::queue<SubCell> stack ({
            SubCell {subcells[i], 1}
        });

        while (not stack.empty()) {
            const auto & subcell = stack.front();
            const auto & local_hexahedron = subcell.hexahedron;

            bool all_local_gauss_points_are_inside = true;
            std::vector<std::pair<LocalCoordinates, Real>> local_gauss_points;
            for (std::size_t gauss_node_id = 0; gauss_node_id < Hexahedron::number_of_gauss_nodes; ++gauss_node_id) {
                const auto &gauss_node = MapVector(Hexahedron::gauss_nodes[gauss_node_id]);
                const auto &gauss_weight = Hexahedron::gauss_weights[gauss_node_id];
                const auto gauss_position = hexa.T(local_hexahedron.T(gauss_node));
                Real fx,fy,fz; // unused
                const auto gauss_cube_id = grid.findCube({gauss_position[0], gauss_position[1], gauss_position[2]}, fx, fy, fz);
                if (gauss_cube_id < 0) {
                    all_local_gauss_points_are_inside = false;
                } else {
                    local_gauss_points.push_back({
                        local_hexahedron.T(gauss_node),
                        gauss_weight*local_hexahedron.jacobian().determinant()
                    });
                }
            }

            if (all_local_gauss_points_are_inside or subcell.subdivision_level == max_number_of_subdivisions) {
                for (const auto & gauss_point : local_gauss_points) {
                    gauss_points[i].push_back(gauss_point);
                }
            } else {
                const auto local_hexahedrons = split_in_local_hexahedrons(subcell.hexahedron);
                for (const auto & local_cell : local_hexahedrons) {
                    stack.push({local_cell, subcell.subdivision_level+1});
                }
            }

            stack.pop();
        }
    }

    return gauss_points;
}

HexahedronElasticForce::HexahedronElasticForce()
: d_youngModulus(initData(&d_youngModulus,
        Real(1000), "youngModulus",
        "Young's modulus of the material",
        true /*displayed_in_GUI*/, false /*read_only_in_GUI*/))
, d_poissonRatio(initData(&d_poissonRatio,
        Real(0.3),  "poissonRatio",
        "Poisson's ratio of the material",
        true /*displayed_in_GUI*/, false /*read_only_in_GUI*/))
, d_number_of_subdivisions(initData(&d_number_of_subdivisions,
        (UNSIGNED_INTEGER_TYPE) 0,  "number_of_subdivisions",
        "Number of subdivisions for the integration method (see 'integration_method').",
        true /*displayed_in_GUI*/, false /*read_only_in_GUI*/))
, d_linear_strain(initData(&d_linear_strain,
        bool(true), "linearStrain",
        "True if the small (linear) strain tensor is used, otherwise the nonlinear Green-Lagrange strain is used.",
        true /*displayed_in_GUI*/, false /*read_only_in_GUI*/))
, d_corotated(initData(&d_corotated,
        bool(true), "corotated",
        "Whether or not to use corotated elements for the strain computation.",
        true /*displayed_in_GUI*/, false /*read_only_in_GUI*/))
, d_ignore_volume_threshold(initData(&d_ignore_volume_threshold,
        Real(0.),  "ignore_volume_threshold",
        "Ignore hexahedrons for which the ratio of the volume inside the boundary against the hexa's volume is "
              "lower than the given threshold",
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
, d_topology_container(initLink(
        "topology_container", "Topology that contains the elements on which this force will be computed."))
, d_integration_grid(initLink(
        "integration_grid", "Sparse grid used for subdivided integration methods (see 'integration_method')."))
{
    d_integration_method.setValue(sofa::helper::OptionsGroup(std::vector<std::string> {
        "Regular", "SubdividedVolume", "SubdividedGauss"
    }));

    sofa::helper::WriteAccessor<Data< sofa::helper::OptionsGroup >> integration_method = d_integration_method;
    integration_method->setSelectedItem((unsigned int) 0);
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

    if (integration_method() == IntegrationMethod::SubdividedVolume or
        integration_method() == IntegrationMethod::SubdividedGauss) {

        if (not d_integration_grid.get()) {
            msg_error() << "Integration method '" << integration_method_as_string()
                        << "' requires a sparse grid topology.";
            set_integration_method(IntegrationMethod::Regular);
        } else {

            if (d_integration_grid.get()->getNbHexas() == 0) {
                msg_error() << "Integration grid '" << d_integration_grid.get()->getPathName()
                            << "' does not contain any hexahedron.";
                set_integration_method(IntegrationMethod::Regular);
            }

            if (d_number_of_subdivisions.getValue() == 0) {
                msg_error() << "Integration method '" << integration_method_as_string()
                            << "' requires the number of subdivisions to be greater than 0";
                set_integration_method(IntegrationMethod::Regular);
            }
        }
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
    const Real volume_ratio_threshold = d_ignore_volume_threshold.getValue();
    std::vector<std::vector<GaussNode>> quadrature_nodes(topology->getNbHexahedra());
    p_quadrature_nodes.resize(topology->getNbHexahedra());
    p_hexahedrons_indices.resize(topology->getNbHexahedra());
    p_ignored_hexahedrons_indices.resize(topology->getNbHexahedra());

    UNSIGNED_INTEGER_TYPE next_hexahedron_id = 0;
    UNSIGNED_INTEGER_TYPE next_ignored_hexahedron_id = 0;
    std::map<float, UNSIGNED_INTEGER_TYPE> volume_ratios;

    for (std::size_t hexa_id = 0; hexa_id < topology->getNbHexahedra(); ++hexa_id) {
        auto   hexa = hexahedron(hexa_id, X);

        // List of pair (Xi, w) where Xi is the local coordinates vector of the gauss point and w is its weight.
        std::vector<std::pair<Vec3, Real>> gauss_points;

        if (integration_method() == IntegrationMethod::Regular) {
            for (std::size_t gauss_node_id = 0; gauss_node_id < Hexahedron::number_of_gauss_nodes; ++gauss_node_id) {
                const auto &gauss_node   = MapVector<3>(Hexahedron::gauss_nodes[gauss_node_id]);
                const auto &gauss_weight = Hexahedron::gauss_weights[gauss_node_id];

                gauss_points.emplace_back(gauss_node, gauss_weight);
            }
        } else {
            auto & grid = *d_integration_grid.get();

            auto gauss_points_per_subcells = recursively_get_subcells_gauss_points(grid, hexa, d_number_of_subdivisions.getValue());
            if (integration_method() == IntegrationMethod::SubdividedVolume) {
                // SubdividedVolume method keeps only the top level hexahedron gauss points, but adjust their weight
                // by integrating the volume of the subcell.
                for (std::size_t gauss_node_id = 0; gauss_node_id < Hexahedron::number_of_gauss_nodes; ++gauss_node_id) {
                    const auto &gauss_node   = MapVector<3>(Hexahedron::gauss_nodes[gauss_node_id]);
                    const auto &gauss_weight = Hexahedron::gauss_weights[gauss_node_id];

                    Real volume = 0.;
                    for (const auto & subcell_gauss_point : gauss_points_per_subcells[gauss_node_id]) {
                        volume += subcell_gauss_point.second;
                    }
                    if (volume > 0) {
                        gauss_points.emplace_back(gauss_node, gauss_weight * volume);
                    }
                }
            } else { /* IntegrationMethod::SubdividedGauss */
                // SubdividedGauss method keeps all gauss points found recursively in the subcells
                for (const auto & gauss_points_in_subcell : gauss_points_per_subcells) {
                    for (const auto & gauss_point : gauss_points_in_subcell) {
                        gauss_points.push_back(gauss_point);
                    }
                }
            }
        }

        // At this point, we have all the gauss points of the hexa and their corrected weight.
        // 1. Ignore hexahedrons for which their internal volume ratio is lower than the given threshold
        bool ignore_hexahedron = false;
        if (integration_method() != IntegrationMethod::Regular) {
            Real full_volume = 0.;
            Real true_volume = 0.;

            for (std::size_t gauss_node_id = 0;
                 gauss_node_id < Hexahedron::number_of_gauss_nodes; ++gauss_node_id) {
                const auto &gauss_node = MapVector<3>(Hexahedron::gauss_nodes[gauss_node_id]);
                const auto &gauss_weight = Hexahedron::gauss_weights[gauss_node_id];
                const auto J = hexa.jacobian(gauss_node);
                const auto detJ = J.determinant();
                full_volume += gauss_weight*detJ;
            }

            for (const auto & gauss_point : gauss_points) {
                const auto J = hexa.jacobian(gauss_point.first);
                const auto &gauss_weight = gauss_point.second;
                const auto detJ = J.determinant();
                true_volume += gauss_weight*detJ;
            }

            const auto volume_ratio = (true_volume / full_volume);
            ignore_hexahedron = (volume_ratio < volume_ratio_threshold);
            if (volume_ratios.find(volume_ratio) == volume_ratios.end())
                volume_ratios[volume_ratio] = 1;
            else
                volume_ratios[volume_ratio] += 1;
        } else if (volume_ratio_threshold > 0) {
            msg_warning() << "Giving a volume threshold value has no effect when using the regular integration method.";
        }

        // 2. We now compute constant values required by the simulation for each of them
        if (not ignore_hexahedron) {
            auto &hexa_quadrature_nodes = p_quadrature_nodes[next_hexahedron_id];
            for (const auto &gauss_point : gauss_points) {
                // Jacobian of the gauss node's transformation mapping from the elementary space to the world space
                const auto J = hexa.jacobian(gauss_point.first);
                const Mat33 Jinv = J.inverse();
                const auto detJ = J.determinant();

                // Derivatives of the shape functions at the gauss node with respect to global coordinates x,y and z
                const Matrix<NumberOfNodes, 3, Eigen::RowMajor> dN_dx = (Jinv.transpose() * Hexahedron::dL(
                    gauss_point.first).transpose()).transpose();

                hexa_quadrature_nodes.push_back(GaussNode({
                                                              gauss_point.second,
                                                              detJ,
                                                              dN_dx,
                                                              Mat33::Identity()
                                                          }));
            }

            p_hexahedrons_indices[next_hexahedron_id] = hexa_id;
            ++next_hexahedron_id;
        } else {
            p_ignored_hexahedrons_indices[next_ignored_hexahedron_id] = hexa_id;
            ++next_ignored_hexahedron_id;
        }
    }

    p_quadrature_nodes.resize(next_hexahedron_id);
    p_hexahedrons_indices.resize(next_hexahedron_id);
    p_ignored_hexahedrons_indices.resize(next_ignored_hexahedron_id);

    if (integration_method() != IntegrationMethod::Regular) {
        for (const auto &p : volume_ratios) {
            msg_info() << p.second << " hexahedrons with a volume ratio of " << p.first;
        }
        msg_info() << p_hexahedrons_indices.size() << " hexahedrons ("
                   << p_ignored_hexahedrons_indices.size() << " ignored on a total of " << topology->getNbHexahedra() << ")";
    }

    // Initialize the stiffness matrix of every hexahedrons
    p_stiffness_matrices.resize(p_hexahedrons_indices.size());
    p_initial_rotation.resize(p_hexahedrons_indices.size(), Mat33::Identity());
    p_current_rotation.resize(p_hexahedrons_indices.size(), Mat33::Identity());

    // Initialize the initial frame of each hexahedron
    if (d_corotated.getValue()) {
        if (not d_linear_strain.getValue()) {
            msg_warning() << "The corotated method won't be computed since nonlinear strain is used.";
        } else {
            for (const auto & hexa_id : p_hexahedrons_indices) {
                Hexahedron hexa = hexahedron(hexa_id, X);
                p_initial_rotation[hexa_id] = hexa.frame();
            }
        }
    }

#ifdef CARIBOU_WITH_OPENMP
    msg_info() << "Using " << omp_get_num_threads() << " threads for computations.";
#endif

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

    static const auto I = Matrix<3,3, Eigen::RowMajor>::Identity();

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
    bool linear = d_linear_strain.getValue();
    recompute_compute_tangent_stiffness = (not linear) and mparams->implicit();
    if (linear) {
        // Small (linear) strain
        sofa::helper::AdvancedTimer::stepBegin("HexahedronElasticForce::addForce");
        #pragma omp parallel for default(none) shared(topology, current_rotation, corotated, x, x0, f)
        for (std::size_t index = 0; index < p_hexahedrons_indices.size(); ++index) {
            const auto & hexa_id = p_hexahedrons_indices[index];
            Hexahedron hexa = hexahedron(hexa_id, x);

            const Mat33 & R0 = initial_rotation[index];
            const Mat33 R0t = R0.transpose();

            Mat33 & R = current_rotation[index];


            // Extract the hexahedron's frame
            if (corotated)
                R = hexa.frame();

            const Mat33 & Rt = R.transpose();

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
            const auto &K = p_stiffness_matrices[index];
            Vec24 F = K * U;

            // Write the forces into the output vector
            i = 0;
            for (const auto &node_id : topology->getHexahedron(hexa_id)) {
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
        sofa::helper::AdvancedTimer::stepEnd("HexahedronElasticForce::addForce");
    } else {
        // Nonlinear Green-Lagrange strain
        sofa::helper::AdvancedTimer::stepBegin("HexahedronElasticForce::addForce");
        const Real youngModulus = d_youngModulus.getValue();
        const Real poissonRatio = d_poissonRatio.getValue();

        const Real l = youngModulus * poissonRatio / ((1 + poissonRatio) * (1 - 2 * poissonRatio));
        const Real m = youngModulus / (2 * (1 + poissonRatio));

    #pragma omp parallel for default(none) shared(topology, current_rotation, corotated, x, x0, f, l, m)
        for (std::size_t index = 0; index < p_hexahedrons_indices.size(); ++index) {
            const auto & hexa_id = p_hexahedrons_indices[index];
            const auto &hexa = topology->getHexahedron(hexa_id);

            Matrix<8, 3, Eigen::RowMajor> U;

            for (size_t i = 0; i < 8; ++i) {
                const auto u = x[hexa[i]] - x0[hexa[i]];
                U(i, 0) = u[0];
                U(i, 1) = u[1];
                U(i, 2) = u[2];
            }

            Matrix<8, 3, Eigen::RowMajor> forces;
            forces.fill(0);
            for (GaussNode &gauss_node : p_quadrature_nodes[index]) {

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
                for (size_t i = 0; i < 8; ++i) {
                    const Vec3 dx = dN_dx.row(i).transpose();
                    const Vec3 f_ = (detJ * w) * F.transpose() * (S * dx);
                    forces(i, 0) += f_[0];
                    forces(i, 1) += f_[1];
                    forces(i, 2) += f_[2];
                }
            }

            for (size_t i = 0; i < 8; ++i) {
                #pragma omp atomic
                f[hexa[i]][0] -= forces.row(i)[0];

                #pragma omp atomic
                f[hexa[i]][1] -= forces.row(i)[1];

                #pragma omp atomic
                f[hexa[i]][2] -= forces.row(i)[2];
            }
        }
        sofa::helper::AdvancedTimer::stepEnd("HexahedronElasticForce::addForce");
    }
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

    if (p_stiffness_matrices.size() != p_hexahedrons_indices.size())
        return;

    if (recompute_compute_tangent_stiffness)
        compute_K();

    auto kFactor = (Real)mparams->kFactorIncludingRayleighDamping(this->rayleighStiffness.getValue());
    sofa::helper::ReadAccessor<Data<VecDeriv>> dx = d_dx;
    sofa::helper::WriteAccessor<Data<VecDeriv>> df = d_df;
    std::vector<Mat33> & current_rotation = p_current_rotation;

    sofa::helper::AdvancedTimer::stepBegin("HexahedronElasticForce::addDForce");
    #pragma omp parallel for default (none) shared (topology, current_rotation, dx, df, kFactor)
    for (std::size_t index = 0; index < p_hexahedrons_indices.size(); ++index) {
        const auto & hexa_id = p_hexahedrons_indices[index];

        const Mat33 & R  = current_rotation[index];
        const Mat33 & Rt = R.transpose();

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
        const auto & K = p_stiffness_matrices[index];
        Vec24 F = K*U*kFactor;

        // Write the forces into the output vector
        i = 0;
        for (const auto & node_id : topology->getHexahedron(hexa_id)) {
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
    sofa::helper::AdvancedTimer::stepEnd("HexahedronElasticForce::addDForce");
}

void HexahedronElasticForce::addKToMatrix(sofa::defaulttype::BaseMatrix * matrix, SReal kFact, unsigned int & /*offset*/)
{
    auto * topology = d_topology_container.get();

    if (!topology)
        return;

    if (recompute_compute_tangent_stiffness)
        compute_K();

    std::vector<Mat33> & current_rotation = p_current_rotation;

    sofa::helper::AdvancedTimer::stepBegin("HexahedronElasticForce::addKToMatrix");

    #pragma omp parallel for default (none) shared (topology, current_rotation, matrix, kFact)
    for (std::size_t index = 0; index < p_hexahedrons_indices.size(); ++index) {
        const auto & hexa_id = p_hexahedrons_indices[index];
        const auto & node_indices = topology->getHexahedron(hexa_id);
        const Mat33 & R  = current_rotation[index];
        const Mat33   Rt = R.transpose();

        const auto & K = p_stiffness_matrices[index];

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

                        #pragma omp critical
                        matrix->add(x, y, k(m,n));
                    }
                }
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

    if (p_stiffness_matrices.size() != p_hexahedrons_indices.size())
        return;

    static const auto I = Matrix<3,3, Eigen::RowMajor>::Identity();
    const Real youngModulus = d_youngModulus.getValue();
    const Real poissonRatio = d_poissonRatio.getValue();

    const Real l = youngModulus * poissonRatio / ((1 + poissonRatio) * (1 - 2 * poissonRatio));
    const Real m = youngModulus / (2 * (1 + poissonRatio));

    sofa::helper::AdvancedTimer::stepBegin("HexahedronElasticForce::compute_k");

    #pragma omp parallel for default (none) shared (topology, l, m)
    for (std::size_t index = 0; index < p_hexahedrons_indices.size(); ++index) {
        auto & K = p_stiffness_matrices[index];
        K.fill(0.);

        for (GaussNode &gauss_node : p_quadrature_nodes[index]) {
            // Jacobian of the gauss node's transformation mapping from the elementary space to the world space
            const Real detJ = gauss_node.jacobian_determinant;

            // Derivatives of the shape functions at the gauss node with respect to global coordinates x,y and z
            const auto dN_dx = gauss_node.dN_dx;

            // Gauss quadrature node weight
            const auto w = gauss_node.weight;

            // Deformation tensor at gauss node
            Mat33 & F = gauss_node.F;

            // Strain tensor at gauss node
            const Mat33 C = F * F.transpose();
            const Mat33 E = 1/2. * (C - I);

            // Stress tensor at gauss node
            const Mat33 S = 2.*m*E + (l * E.trace() * I);

            // Computation of the tangent-stiffness matrix
            for (std::size_t i = 0; i < 8; ++i) {
                // Derivatives of the ith shape function at the gauss node with respect to global coordinates x,y and z
                const Vec3 dxi = dN_dx.row(i).transpose();

                for (std::size_t j = 0; j < 8; ++j) {
                    // Derivatives of the jth shape function at the gauss node with respect to global coordinates x,y and z
                    const Vec3 dxj = dN_dx.row(j).transpose();

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
    sofa::helper::AdvancedTimer::stepEnd("HexahedronElasticForce::compute_k");
}

const Eigen::SparseMatrix<HexahedronElasticForce::Real> & HexahedronElasticForce::K() {
    if (not K_is_up_to_date) {
        const sofa::helper::ReadAccessor<Data<VecCoord>> X = this->mstate->readRestPositions();
        const auto nDofs = X.size() * 3;
        p_K.resize(nDofs, nDofs);
        p_K.setZero();
        p_K.reserve(Eigen::VectorXi::Constant(nDofs, 24));

        auto *topology = d_topology_container.get();

        if (topology) {

            const std::vector<Mat33> &current_rotation = p_current_rotation;

            for (std::size_t index = 0; index < p_hexahedrons_indices.size(); ++index) {
                const auto & hexa_id = p_hexahedrons_indices[index];
                const auto &node_indices = topology->getHexahedron(hexa_id);
                const Mat33 &R = current_rotation[index];
                const Mat33 Rt = R.transpose();

                const auto &Ke = p_stiffness_matrices[index];

                for (size_t i = 0; i < 8; ++i) {
                    for (size_t j = 0; j < 8; ++j) {
                        Mat33 k = Ke.block(i, j, 3, 3);

                        k = -1. * R * k * Rt;

                        for (unsigned char m = 0; m < 3; ++m) {
                            for (unsigned char n = 0; n < 3; ++n) {
                                const auto x = node_indices[i] * 3 + m;
                                const auto y = node_indices[j] * 3 + n;

                                const auto v = p_K.coeff(x, y);
                                p_K.coeffRef(x, y) = v + k(m, n);
                            }
                        }
                    }
                }
            }
        }
        p_K.makeCompressed();
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

void HexahedronElasticForce::computeBBox(const sofa::core::ExecParams* params, bool onlyVisible)
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
    for (const auto & hexa_id : p_hexahedrons_indices) {
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
    for (const auto & hexa_id : p_ignored_hexahedrons_indices) {
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

} // namespace SofaCaribou::GraphComponents::forcefield
