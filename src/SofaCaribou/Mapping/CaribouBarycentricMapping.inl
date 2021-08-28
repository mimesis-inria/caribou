#include <SofaCaribou/Mapping/CaribouBarycentricMapping.h>

DISABLE_ALL_WARNINGS_BEGIN
#include <sofa/version.h>
#include <sofa/core/Mapping.inl>
#include <sofa/simulation/Node.h>
#include <sofa/core/visual/VisualParams.h>
DISABLE_ALL_WARNINGS_END

#if (defined(SOFA_VERSION) && SOFA_VERSION < 201200)
namespace sofa { using Index = unsigned int; }
#endif

#if (defined(SOFA_VERSION) && SOFA_VERSION < 210600)
namespace sofa::type {
using RGBAColor = ::sofa::helper::types::RGBAColor;
using Vector3 = ::sofa::defaulttype::Vector3;
}
#endif

namespace SofaCaribou::mapping {

template<typename Element, typename MappedDataTypes>
CaribouBarycentricMapping<Element, MappedDataTypes>::CaribouBarycentricMapping()
: d_topology(initLink("topology", "Topology that contains the embedding (parent) elements."))
{
}

template<typename Element, typename MappedDataTypes>
void CaribouBarycentricMapping<Element, MappedDataTypes>::init() {
    using CaribouTopology = SofaCaribou::topology::CaribouTopology<Element>;

    // Get the topology if it is not already provided by the user
    if (d_topology.empty()) {
        auto * context_node = dynamic_cast<sofa::simulation::Node*>(this->getContext());
        auto * parent_context_node = (context_node->getNbParents() > 0) ? context_node->getFirstParent() : nullptr;
        CaribouTopology * topology = nullptr;

        // Search first in the parent context node (if any)
        if (parent_context_node) {
            topology = parent_context_node->getContext()->get<CaribouTopology>();
        }

        // If not found in the parent context, search in the current context
        if (not topology and (context_node != parent_context_node)) {
            topology = context_node->getContext()->get<CaribouTopology>();
        }

        if (not topology) {
            msg_error()
                    << "Could not find a topology container in the parent context node, or the current context node. "
                    << "Please add a '" << CaribouTopology::GetClass()->typeName << "' component to your scene, and set its "
                    << "complete path in the '" << d_topology.getName() << "' data attribute.";
            return;
        }
        // We got one, use it
        d_topology.set(topology);
        msg_info() << "Topology automatically found at '" << topology->getPathName() << "'";
    }

    // Quick sanity checks
    if (not d_topology->domain()) {
        msg_error() << "The topology '" << d_topology->getPathName() << "' doesn't seem to have been initialized, or "
                    << "might not contain any elements.";
        return;
    }

    if (this->getFromModel()->getSize() == 0) {
        msg_error() << "The parent mechanical object '" << this->getFromModel()->getPathName() << "' doesn't contain any nodes.";
        return;
    }

    if (this->getToModel()->getSize() == 0) {
        msg_error() << "The mapped mechanical object '" << this->getToModel()->getPathName() << "' doesn't contain any nodes.";
        return;
    }

    // Initialize the rest positions
    if (this->getToModel()->readRestPositions().empty()) {
        auto mapped_rest_positions = this->getToModel()->writeOnlyRestPositions();
        auto mapped_positions = this->getToModel()->readPositions();
        mapped_rest_positions.resize(mapped_positions.size());
        for (sofa::Index i = 0; i < mapped_positions.size(); ++i) {
            mapped_rest_positions[i] = mapped_positions[i];
        }
    }

    // Create the barycentric container
    auto mapped_rest_positions =
        Eigen::Map<const Eigen::Matrix<Scalar, Eigen::Dynamic, Dimension, Eigen::RowMajor>> (
            this->getToModel()->readRestPositions()[0].ptr(),
            this->getToModel()->readRestPositions().size(),
            Dimension
        );
    p_barycentric_container.reset(new caribou::topology::BarycentricContainer<Domain>(d_topology->domain()->embed(mapped_rest_positions)));

    if (not p_barycentric_container->outside_nodes().empty()) {
        const auto n = p_barycentric_container->outside_nodes().size();
        msg_error() << n << " / " << mapped_rest_positions.rows() << " "
                      << "mapped nodes were found outside of the embedding topology, i.e., they are not "
                      << "contained inside any elements. Every nodes must be found inside the domain.";
        return;
    }

    // Precompute the mapping matrix (derivative of the mapping function w.r.t rest position).
    p_J.setZero(); // Let's clear it just in case init was called previously
    p_J.resize(this->getToModel()->getSize(), this->getFromModel()->getSize());
    std::vector<Eigen::Triplet<Scalar>> entries;
    entries.reserve(this->getFromModel()->getSize());
    const auto * domain = d_topology->domain();
    const auto & barycentric_points = p_barycentric_container->barycentric_points();
    for (std::size_t i = 0; i < barycentric_points.size(); ++i) {
        const auto & bp = barycentric_points[i];
        const auto & node_indices = domain->element_indices(bp.element_index);
        const auto e = domain->element(bp.element_index);
        const auto L = e.L(bp.local_coordinates); // Shape functions at barycentric coordinates
        for (Eigen::Index j = 0; j < node_indices.rows(); ++j) {
            const auto & node_index = node_indices[j];
            entries.emplace_back(i, node_index, L[j]);
        }
    }
    p_J.setFromTriplets(entries.begin(), entries.end());

    // This is needed to map the initial nodal positions and velocities (and, optionally, rest positions).
    Inherit1 ::init();
}

template<typename Element, typename MappedDataTypes>
void CaribouBarycentricMapping<Element, MappedDataTypes>::apply(const sofa::core::MechanicalParams * /*mparams*/,
                                                                MappedDataVecCoord & data_output_mapped_position,
                                                                const DataVecCoord & data_input_position) {
    // Sanity check
    if (not p_barycentric_container or not p_barycentric_container->outside_nodes().empty()) {
        return;
    }

    // Container (parent) nodes
    auto input_position = sofa::helper::ReadAccessor<DataVecCoord>(data_input_position);
    auto positions =
            Eigen::Map<const Eigen::Matrix<Scalar, Eigen::Dynamic, Dimension, Eigen::RowMajor>> (
                    input_position[0].ptr(),
                    input_position.size(),
                    Dimension
            );

    // Mapped (embedded) nodes
    auto output_mapped_position = sofa::helper::WriteOnlyAccessor<MappedDataVecCoord>(data_output_mapped_position);
    auto mapped_positions =
            Eigen::Map<Eigen::Matrix<MappedScalar, Eigen::Dynamic, MappedDimension, Eigen::RowMajor>> (
                    output_mapped_position[0].ptr(),
                    output_mapped_position.size(),
                    MappedDimension
            );

    // Interpolate the parent positions onto the mapped nodes
    p_barycentric_container->interpolate(positions, mapped_positions);
}

template<typename Element, typename MappedDataTypes>
void CaribouBarycentricMapping<Element, MappedDataTypes>::applyJ(const sofa::core::MechanicalParams * /*mparams*/,
                                                                 MappedDataVecDeriv & data_output_mapped_velocity,
                                                                 const DataVecDeriv & data_input_velocity) {
    // Sanity check
    if (not p_barycentric_container or not p_barycentric_container->outside_nodes().empty()) {
        return;
    }

    // Container (parent) nodes
    auto input_velocities = sofa::helper::ReadAccessor<DataVecDeriv>(data_input_velocity);
    auto velocities =
            Eigen::Map<const Eigen::Matrix<Scalar, Eigen::Dynamic, Dimension, Eigen::RowMajor>> (
                    input_velocities[0].ptr(),
                    input_velocities.size(),
                    Dimension
            );

    // Mapped (embedded) nodes
    auto output_mapped_position = sofa::helper::WriteOnlyAccessor<MappedDataVecDeriv>(data_output_mapped_velocity);
    auto mapped_velocities =
            Eigen::Map<Eigen::Matrix<MappedScalar, Eigen::Dynamic, MappedDimension, Eigen::RowMajor>> (
                    output_mapped_position[0].ptr(),
                    output_mapped_position.size(),
                    MappedDimension
            );

    // Interpolate the parent positions onto the mapped nodes
    p_barycentric_container->interpolate(velocities, mapped_velocities);
}

template<typename Element, typename MappedDataTypes>
void CaribouBarycentricMapping<Element, MappedDataTypes>::applyJT(const sofa::core::MechanicalParams * /*mparams*/,
                                                                  DataVecDeriv & data_output_force,
                                                                  const MappedDataVecDeriv & data_input_mapped_force) {
    // Sanity check
    if (not p_barycentric_container or not p_barycentric_container->outside_nodes().empty()) {
        return;
    }

    // Container (parent) nodes
    auto output_forces = sofa::helper::WriteAccessor<DataVecDeriv>(data_output_force);
    auto forces =
            Eigen::Map<Eigen::Matrix<Scalar, Eigen::Dynamic, Dimension, Eigen::RowMajor>> (
                    output_forces[0].ptr(),
                    output_forces.size(),
                    Dimension
            );

    // Mapped (embedded) nodes
    auto input_mapped_forces = sofa::helper::ReadAccessor<MappedDataVecDeriv>(data_input_mapped_force);
    auto mapped_forces =
            Eigen::Map<const Eigen::Matrix<MappedScalar, Eigen::Dynamic, MappedDimension, Eigen::RowMajor>> (
                    input_mapped_forces[0].ptr(),
                    input_mapped_forces.size(),
                    MappedDimension
            );

    // Inverse mapping using the transposed of the Jacobian matrix
    forces.noalias() += p_J.transpose() * mapped_forces;

}

template<typename Element, typename MappedDataTypes>
void CaribouBarycentricMapping<Element, MappedDataTypes>::applyJT(const sofa::core::ConstraintParams * /*cparams*/,
                                                                  CaribouBarycentricMapping::DataMapMapSparseMatrix & data_output_jacobian,
                                                                  const CaribouBarycentricMapping::MappedDataMapMapSparseMatrix & data_input_mapped_jacobian) {
    using SparseMatrix = decltype(p_J);
    using StorageIndex = typename SparseMatrix::StorageIndex;
    static constexpr auto J_is_row_major = SparseMatrix::IsRowMajor;

    // Sanity check
    if (not p_barycentric_container or not p_barycentric_container->outside_nodes().empty()) {
        return;
    }

    // With the mapping matrix J, every columns j of a given row i represent the jth shape
    // value evaluated at the ith mapped node. We want to compute H = Jt H_m where H_m is
    // the constraint matrix of the mapped DOFs. Both Jt and H_m are sparse matrices.
    // However, H and H_m are a special kind of sparse matrices having a map of map of "Deriv"
    // vector.
    auto output_jacobian = sofa::helper::WriteAccessor<DataMapMapSparseMatrix>(data_output_jacobian);
    auto input_mapped_jacobian = sofa::helper::ReadAccessor<MappedDataMapMapSparseMatrix>(data_input_mapped_jacobian);

    const auto rowEnd = input_mapped_jacobian->end();
    for (auto rowIt = input_mapped_jacobian->begin(); rowIt != rowEnd; ++rowIt) {
        auto colItEnd = rowIt.end();
        auto colIt = rowIt.begin();
        if (colIt == colItEnd) continue;

        auto output_row = output_jacobian->writeLine(rowIt.index());
        for ( ; colIt != colItEnd; ++colIt) {
            auto mapped_node_index  =  colIt.index();
            auto mapped_value = colIt.val();

            if constexpr (J_is_row_major) {
                for (typename SparseMatrix::InnerIterator it(p_J, mapped_node_index); it; ++it) {
                    const auto node_index = it.col();
                    const auto shape_value = it.value();
                    output_row.addCol(node_index, mapped_value*shape_value);
                }
            } else {
                // J is column major
                for (int node_index=0; node_index<p_J.outerSize(); ++node_index) {
                    const auto & start = p_J.outerIndexPtr()[node_index];
                    const auto & end   = p_J.outerIndexPtr()[node_index+1];
                    const auto id = p_J.data().searchLowerIndex(start, end, static_cast<StorageIndex>(mapped_node_index));
                    if ((id<end) && (p_J.data().index(id)==static_cast<StorageIndex>(mapped_node_index))) {
                        const auto shape_value = p_J.data().value(id);
                        output_row.addCol(node_index, mapped_value*shape_value);
                    }
                }
            }
        }
    }
}

template<typename Element, typename MappedDataTypes>
void CaribouBarycentricMapping<Element, MappedDataTypes>::draw(const sofa::core::visual::VisualParams *vparams) {
    if ( !vparams->displayFlags().getShowMappings() ) return;

    if (not p_barycentric_container) return;

    using Color = sofa::type::RGBAColor;

    // Draw outside points
    if (not p_barycentric_container->outside_nodes().empty()) {
        const auto mapped_rest_positions = this->getToModel()->readRestPositions();
        const auto & outside_nodes = p_barycentric_container->outside_nodes();
        std::vector<sofa::type::Vector3> outside_nodes_positions (outside_nodes.size());
        for (std::size_t i = 0; i < outside_nodes.size(); ++i) {
            const auto & node_index = outside_nodes[i];
            for (std::size_t j = 0; j < MappedDimension; ++j) {
                outside_nodes_positions[i][j] =    mapped_rest_positions[node_index][j];
            }
        }
        vparams->drawTool()->drawPoints(outside_nodes_positions, 5, Color::red());
    }
}

template<typename Element, typename MappedDataTypes>
template<typename Derived>
auto CaribouBarycentricMapping<Element, MappedDataTypes>::canCreate(Derived * o,
                                                                    sofa::core::objectmodel::BaseContext * context,
                                                                    sofa::core::objectmodel::BaseObjectDescription * arg) -> bool {
    using namespace sofa::core::objectmodel;
    using CaribouTopology = SofaCaribou::topology::CaribouTopology<Element>;

    std::string requested_element_type = arg->getAttribute( "template", "");
    std::string this_element_type = Derived::templateName(o);

    // to lower
    std::string requested_element_type_lower = requested_element_type, this_element_type_lower = this_element_type;
    std::transform(requested_element_type.begin(), requested_element_type.end(), requested_element_type_lower.begin(), [](unsigned char c){ return std::tolower(c); });
    std::transform(this_element_type.begin(), this_element_type.end(), this_element_type_lower.begin(), [](unsigned char c){ return std::tolower(c); });

    // Simplest case, the user set the template argument and it is equal to this template specialization
    if (requested_element_type_lower == this_element_type_lower) {
        return Inherit1::canCreate(o, context, arg);
    }

    // The user set the template argument, but it doesn't match this template specialization
    if (not requested_element_type.empty()) {
        arg->logError("Requested element type is not '"+this_element_type+"'.");
        return false;
    }
    // The user didn't set any template argument. Let's try to deduce it.

    // The user specified a path to a CaribouTopology, let's check if the Element type matches
    std::string topology_path = arg->getAttribute("topology", "");
    if (not topology_path.empty()) {
        topology_path = topology_path.substr(1); // removes the "@"
        auto topology = context->getRootContext()->get<CaribouTopology>(topology_path);
        if (not topology) {
            arg->logError("The element type of the given topology '" + topology_path + "' is not '"+this_element_type+"'.");
            return false;
        }
        // The topology element type matched this mapping, all good here
        return Inherit1::canCreate(o, context, arg);
    }

    // The user hasn't specified any path to a topology, nor he has manually set the Element type using the template argument
    // Let's try to find a matching topology in the parent or current context.
    auto * context_node = dynamic_cast<sofa::simulation::Node*>(context);
    auto * parent_context_node = (context_node->getNbParents() > 0) ? context_node->getFirstParent() : nullptr;
    CaribouTopology * topology = nullptr;

    // Search first in the parent context node (if any)
    if (parent_context_node) {
        topology = parent_context_node->getContext()->get<CaribouTopology>();
    }

    // If not found in the parent context, search in the current context
    if (not topology and (context_node != parent_context_node)) {
        topology = context_node->getContext()->get<CaribouTopology>();
    }

    if (not topology) {
        arg->logError("Cannot deduce the element type from any topology in the parent or current context nodes.");
        return false;
    }

    // We found a topology with a matching element type, let's try it. The user will be noticed of this automatic choice
    // made for him in the init method
    return Inherit1::canCreate(o, context, arg);
}

} // namespace SofaCaribou::mapping

