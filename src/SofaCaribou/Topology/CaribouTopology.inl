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
: d_position(initData(&d_position,
    "position",
    "Position vector of the domain's nodes."))
, d_indices(initData(&d_indices,
    "indices",
    "Node indices (w.r.t the position vector) of each elements."))
{}

template <typename Element>
void CaribouTopology<Element>::attachDomain(const caribou::topology::Domain<Element, PointID> * domain) {
    this->p_domain = domain;

    using namespace sofa::helper;
    auto indices = WriteOnlyAccessor<Data<sofa::type::vector<sofa::type::fixed_array<PointID, NumberOfNodes>>>>(d_indices);

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
    auto indices = ReadAccessor<Data<sofa::type::vector<sofa::type::fixed_array<PointID, NumberOfNodes>>>>(d_indices);
    if (indices.empty()) {
        msg_warning() << "Initializing the topology from an empty set of indices. Make sure you fill the "
                      << "'" << d_indices.getName() << "' data attribute with a vector of node indices.";
        return;
    }

    auto positions = ReadAccessor<Data<VecCoord>>(d_position);
    if (positions.empty()) {
        msg_warning() << "Initializing the topology from a set of indices, but without any node positions "
                      << "vector. Make sure you fill the " << "'" << d_position.getName() << "' data attribute "
                      << "with a vector of node positions.";
        return;
    }

    std::size_t invalid_elements = 0;
    for (std::size_t element_id = 0; element_id < indices.size(); ++element_id) {
        for (const auto & node_index : indices[element_id]) {
            if (node_index >= positions.size()) {
                ++invalid_elements;
                break;
            }
        }
    }
    if (invalid_elements > 0) {
        msg_warning() << "There are " << invalid_elements << " elements containing one or more node index greater than "
                      << "the size of the position vector. Make sure the node indices are relative to the position "
                      << "vector set in the data attribute " << "'" << d_position.getName() << "'.";
        return;
    }

    // Get the mechanical object's position vector and create the internal Mesh instance

    const auto number_of_nodes = positions.size();
    std::vector<typename caribou::topology::Mesh<Dimension>::WorldCoordinates> mesh_positions_vector;
    mesh_positions_vector.reserve(number_of_nodes);
    for (std::size_t i = 0; i < number_of_nodes; ++i) {
        const auto & n = positions[i];
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
        if (not d_indices.getValue().empty() and not d_position.isSet()) {
            auto * context = this->getContext();
            auto state = context->template get<sofa::core::State<DataTypes>>(BaseContext::Local);
            if (state and state->findData("rest_position")) {
                d_position.setParent(state->findData("rest_position"));
                msg_info() << "Automatically found the nodes positions from'" << state->findData("rest_position")->getLinkPath() << "'.";
            }
        }

        this->initializeFromIndices();
    }


}

} // namespace SofaCaribou::topology
