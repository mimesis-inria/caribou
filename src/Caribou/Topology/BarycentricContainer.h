#pragma once

#include <memory>
#include <vector>
#include <unordered_map>
#include <map>

#include <Caribou/config.h>
#include <Caribou/macros.h>
#include <Caribou/Constants.h>
#include <Caribou/Topology/Domain.h>
#include <Caribou/Topology/HashGrid.h>

namespace caribou::topology {

/**
 * \class BarycentricContainer
 *
 * The BarycentricContainer class allows to embed one or more meshes into a container mesh domain.
 * The nodes of embedded meshes do not need to match the nodes of the container domain. However, their nodes must
 * be found in exactly one of the element of the container, or at the boundary between elements (for example, a
 * node of an embedded mesh lying on a face between two elements of the container is valid). If a node of one of
 * the embedded meshes is lying outside of the container domain (ie is not located inside any of the container
 * elements), the BarycentricContainer will ignore it. The list of ignored nodes can be retrieve once an embedded
 * mesh has been added.
 *
 * @tparam ContainerElement The type of element of the container domain.
 * @tparam NodeIndex The type of integer used for a node index of the container domain.
 */
 template <typename ContainerElement, typename NodeIndex>
class BarycentricContainer {
public:

    static constexpr INTEGER_TYPE Dimension = ContainerElement::Dimension;
    using ElementIndex = INTEGER_TYPE;
    using LocalCoordinates = typename ContainerElement::LocalCoordinates;
    using WorldCoordinates = typename ContainerElement::WorldCoordinates;
    using HashGridT = HashGrid<ContainerElement>;

    /**
     * A barycentric point is a structure that contains the element index and
     * the local coordinates within that element of the point in the world space.
     *
     * \note If the point does not exists in any of the containing elements, the index will be -1.
     */
    struct BarycentricPoint {
        ElementIndex element_index;
        LocalCoordinates local_coordinates;
    };

    /** Default constructor is not permitted */
    BarycentricContainer() = delete;

    /**
     * Construct the container from the given domain. This will create an HashGrid class instance to be able
     * to quickly retrieve the container element of a given world position. The size of the HashGrid cells will
     * be set to the mean size of the container elements.
     * @param container_domain The mesh domain that will contain the embedded meshes.
     */
    explicit BarycentricContainer(const Domain<ContainerElement, NodeIndex> * container_domain)
    : p_container_domain(container_domain) {
        if (container_domain->number_of_elements() == 0) {
            throw std::runtime_error("Trying to create a barycentric container from an empty domain.");
        }

        // Get the mean size of the elements
        WorldCoordinates H_mean = WorldCoordinates::Zero();
        for (UNSIGNED_INTEGER_TYPE element_id = 0; element_id < container_domain->number_of_elements(); ++element_id) {
            const ContainerElement element = container_domain->element(element_id);
            const auto element_nodes = element.nodes();
            const auto min = element_nodes.colwise().minCoeff();
            const auto max = element_nodes.colwise().maxCoeff();
            H_mean += (max - min);
        }
        H_mean /= container_domain->number_of_elements();

        // Create the Hash grid
        p_hash_grid = std::make_unique<HashGridT>(H_mean.maxCoeff(), container_domain->number_of_elements());

        for (UNSIGNED_INTEGER_TYPE element_id = 0; element_id < container_domain->number_of_elements(); ++element_id) {
            p_hash_grid->add(container_domain->element(element_id), element_id);
        }
    }

    /**
     * Get an element that contains the given point (in world coordinates) and its local coordinates within this element.
     *
     * @param p The queried point in world coordinates.
     * @return A pair where the first member is an element index within the containing domain, or -1 if the given point
     *         was not inside any elements. The second member returned by this method is the local coordinates of the
     *         point within the containing element. See BarycentricContainer::BarycentricPoint for more details.
     *
     * \warning If the given point lies directly between two or more elements (for example, if the point is on a face,
     *          an edge or a node of the domain), the first element found containing the point will be return.
     */
    auto barycentric_point(const WorldCoordinates & p) const -> BarycentricPoint {
        if (not p_hash_grid) {
            return {-1, LocalCoordinates::Zero()};
        }

        // List of candidate elements that could contain the point
        const auto candidate_element_indices = p_hash_grid->get(p);

        // List of elements that do contain the point
        std::vector<std::pair<INTEGER_TYPE, LocalCoordinates>> elements;
        elements.reserve(candidate_element_indices.size());

        for (const auto & element_index : candidate_element_indices) {
            const ContainerElement e = p_container_domain->element(element_index);
            const LocalCoordinates local_coordinates = e.local_coordinates(p);
            if (e.contains_local(local_coordinates)) {
                // Found one, let's return it
                return {static_cast<ElementIndex>(element_index), local_coordinates};
            }
        }

        // No elements are containing the queried point
        return {-1, LocalCoordinates::Zero()};
    }

    /**
     * Get an element that contains the given embedded point and its local coordinates within this element.
     * @param embedded_mesh The embedded mesh in which the queried point belongs
     * @param embedded_node_index The node indices of the queried point inside the embedded mesh.
     * @return A pair where the first member is an element index within the containing domain, or -1 if the given point
     *         was not inside any elements. The second member returned by this method is the local coordinates of the
     *         point within the containing element. See BarycentricContainer::BarycentricPoint for more details.
     *
     * \warning If the given point lies directly between two or more elements (for example, if the point is on a face,
     *          an edge or a node of the domain), the first element found containing the point will be return.
     */
    auto barycentric_point(const Mesh<Dimension> & embedded_mesh, const NodeIndex & embedded_node_index) const -> BarycentricPoint {
        const auto mesh_entry = p_embedded_meshes.find(&embedded_mesh);
        if (mesh_entry == p_embedded_meshes.end()) {
            // The embedded mesh is not registered
            return {-1, LocalCoordinates::Zero()};
        }

        const auto node_entry = mesh_entry->second.find(embedded_node_index);
        if (node_entry == mesh_entry->second.end()) {
            // The queried point is lying outside the containing mesh
            return {-1, LocalCoordinates::Zero()};
        }

        // Return the barycentric point registered for the queried point
        return node_entry->second;
    }

    /**
     * Add an embedded mesh to the container.
     * @param mesh The embedded mesh.
     * @return The list of node indices of the embedded mesh for which their positions
     *         was found outside of the container domain. If the embedded mesh lies
     *         completely inside the containing domain, this will be empty.
     */
    auto add_embedded_mesh(const Mesh<Dimension> * embedded_mesh) -> std::vector<UNSIGNED_INTEGER_TYPE> {
        std::vector<UNSIGNED_INTEGER_TYPE> outside_nodes;
        outside_nodes.reserve(std::floor(embedded_mesh->number_of_nodes() / 10)); // Reserve 10% of the mesh size for outside nodes

        std::map<UNSIGNED_INTEGER_TYPE, BarycentricPoint> inside_nodes;

        for (UNSIGNED_INTEGER_TYPE node_id = 0; node_id < embedded_mesh->number_of_nodes(); ++node_id) {
            BarycentricPoint b_point = barycentric_point(embedded_mesh->position(node_id));
            if (b_point.element_index < 0) {
                outside_nodes.emplace_back(node_id);
            } else {
                inside_nodes.emplace_hint(inside_nodes.end(), node_id, b_point);
            }
        }

        p_embedded_meshes[embedded_mesh] = inside_nodes;
        return outside_nodes;
    }
private:
    const Domain<ContainerElement, NodeIndex> * p_container_domain;
    std::unordered_map<const Mesh<Dimension> *, std::map<UNSIGNED_INTEGER_TYPE, BarycentricPoint>> p_embedded_meshes;
    std::unique_ptr<HashGridT> p_hash_grid;
};


 // Template deduction guides
 template <typename ContainerElement, typename NodeIndex>
 BarycentricContainer(const Domain<ContainerElement, NodeIndex> *) -> BarycentricContainer<ContainerElement, NodeIndex>;

} // namespace caribou::topology