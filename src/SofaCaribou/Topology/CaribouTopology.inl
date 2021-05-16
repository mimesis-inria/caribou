#pragma once

#include <SofaCaribou/Topology/CaribouTopology.h>

DISABLE_ALL_WARNINGS_BEGIN
#include <sofa/version.h>
DISABLE_ALL_WARNINGS_END

#if (defined(SOFA_VERSION) && SOFA_VERSION < 201299)
namespace sofa { using Index = unsigned int; }
#endif

namespace SofaCaribou::topology {

template <typename Element>
CaribouTopology<Element>::CaribouTopology ()
: d_state (initLink(
    "state",
    "Mechanical state (object) containing the position of each nodes of this topology."))
, d_indices(initData(&d_indices,
    "indices",
    "Contains for each element, a list of each element node indices relative to the mechanical state position vector."))
{}

template <typename Element>
void CaribouTopology<Element>::attachDomain(const caribou::topology::Domain<Element, PointID> * domain) {
    this->p_domain = domain;

    using namespace sofa::helper;
    auto indices = WriteOnlyAccessor<Data<vector<fixed_array<PointID, NumberOfNodes>>>>(d_indices);

    const auto number_of_elements = domain->number_of_elements();
    indices.resize(number_of_elements);

    for (sofa::Index element_id = 0; element_id < number_of_elements; ++element_id) {
        const auto element_indices = domain->element_indices(element_id);
        for (sofa::Index node_id = 0; node_id < NumberOfNodes; ++node_id) {
            indices[element_id][node_id] = static_cast<PointID>(element_indices[node_id]);
        }
    }
}

template <typename Element>
void CaribouTopology<Element>::initializeFromIndices() {
    using namespace sofa::helper;

    // Sanity checks
    auto indices = ReadAccessor<Data<vector<fixed_array<PointID, NumberOfNodes>>>>(d_indices);
    if (indices.empty()) {
        msg_warning() << "Initializing the topology from an empty set of indices. Make sure you fill the "
                      << "'" << d_indices.getName() << "' data attribute with a vector of node indices.";
        return;
    }

    if (not this->d_state.get()) {
        msg_error() << "Initializing the topology from a set of indices requires a mechanical state. The data attribute "
                    << "'" << d_state.getName() << "' can be used to link this component to a mechanical object.";
        return;
    }

    // Get the mechanical object's position vector and create the internal Mesh instance

    const auto & mo_positions = this->d_state->readRestPositions();
    const auto number_of_nodes = mo_positions.size();
    std::vector<typename caribou::topology::Mesh<Dimension>::WorldCoordinates> mesh_positions_vector;
    mesh_positions_vector.reserve(number_of_nodes);
    for (std::size_t i = 0; i < number_of_nodes; ++i) {
        const auto & n = mo_positions[i];
        if constexpr (Dimension == 1) {
            mesh_positions_vector.emplace_back(n[0]);
        } else if constexpr (Dimension == 2) {
            mesh_positions_vector.emplace_back(n[0], n[1]);
        } else {
            mesh_positions_vector.emplace_back(n[0], n[1], n[2]);
        }
    }
    this->p_mesh.reset(new caribou::topology::Mesh<Dimension>(mesh_positions_vector));


    // Read the indices and create the internal Domain instance
    const auto * indices_ptr = indices[0].data();
    this->p_domain = this->p_mesh->template add_domain<Element, PointID>(indices_ptr, indices.size(), NumberOfNodes);
}

template<typename Element>
void CaribouTopology<Element>::init() {
    using namespace sofa::core::objectmodel;

    if (d_indices.isSet()) {
        if (p_domain != nullptr) {
            msg_warning()
                    << "Initializing the topology from a set of indices, while a domain has already been attached. The "
                       "latter will be overridden.";
        }

        // If some indices are set, but no mechanical state is linked, let's try to find one in the current context node
        if (not d_indices.getValue().empty() and d_state.get() == nullptr) {
            auto * context = this->getContext();
            auto state = context->template get<sofa::core::State<DataTypes>>(BaseContext::Local);
            if (state) {
                d_state.set(state);
                msg_info() << "Automatically found the mechanical state '" << d_state.get()->getPathName() << "'.";
            }
        }

        this->initializeFromIndices();
    }


}

} // namespace SofaCaribou::topology