#pragma once

#include <memory>
#include <vector>
#include <unordered_map>

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
 template <typename Domain>
class BarycentricContainer {
public:
    using Mesh = typename Domain::MeshType;
    using ContainerElement = typename Domain::ElementType;

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
    explicit BarycentricContainer(const Domain * container_domain)
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
    auto barycentric_point(const BaseMesh & embedded_mesh, const UNSIGNED_INTEGER_TYPE & embedded_node_index) const -> BarycentricPoint {
        const auto mesh_entry = p_embedded_meshes.find(&embedded_mesh);
        if (mesh_entry == p_embedded_meshes.end()) {
            // The embedded mesh is not registered
            throw std::runtime_error("The embedded mesh hasn't been registered to the barycentric container. "
                                     "Make sure you called the add_embedded_mesh method before calling this one.");
        }

        const std::vector<BarycentricPoint> & barycentric_points = mesh_entry->second;
        caribou_assert(embedded_node_index < barycentric_points.size());

        // Return the barycentric point registered for the queried point
        return barycentric_points[embedded_node_index];
    }

    /**
     * Get the barycentric points of a set of positions embedded inside the container mesh.
     * @tparam Derived NXD Eigen matrix representing the D dimensional coordinates of the N embedded positions.
     * @param embedded_points The positions (in world coordinates) embedded in the container mesh for which the barycentric
     *                        points have to be found.
     * @return Barycentric points (element_id, local coordinates) of every embedded positions passed as argument.
     * \note When an embedded position lie completly outside the container mesh (i.e. it is not contained inside any
     *       container element), the element id of the barycentric point for that position will be -1.
     */
    template <typename Derived>
    auto barycentric_points(const Eigen::MatrixBase<Derived> & embedded_points) const -> std::vector<BarycentricPoint> {
        static_assert(Eigen::MatrixBase<Derived>::ColsAtCompileTime == Dimension or Eigen::MatrixBase<Derived>::ColsAtCompileTime == Eigen::Dynamic);

        if (Eigen::MatrixBase<Derived>::ColsAtCompileTime == Eigen::Dynamic and embedded_points.cols() != Dimension) {
            throw std::runtime_error("Trying to get the barycentric coordinates of " +
                                     std::to_string(embedded_points.cols()) + "D points from a " +
                                     std::to_string(Dimension) + "D mesh");
        }

        std::vector<BarycentricPoint> barycentric_points;
        barycentric_points.reserve(embedded_points.rows());

        const auto number_of_embedded_points = embedded_points.rows();
        for (Eigen::Index node_id = 0; node_id < number_of_embedded_points; ++node_id) {
            barycentric_points.emplace_back(barycentric_point(embedded_points.row(node_id).template cast<typename WorldCoordinates::Scalar>()));
        }

        return barycentric_points;
    }

    /**
     * Add an embedded mesh to the container. If one or more node of the embedded mesh are found outside of the
     * containing domain (ie, no elements containing the node can be found), the operation will fail and the
     * list of nodes found outside the domain are returned.
     *
     * @param mesh The embedded mesh.
     * @return The list of node indices of the embedded mesh for which their positions
     *         was found outside of the container domain. If the embedded mesh lies
     *         completely inside the containing domain, this will be empty.
     */
    template <typename EmbeddedMesh>
    auto add_embedded_mesh(const EmbeddedMesh * embedded_mesh) -> std::vector<UNSIGNED_INTEGER_TYPE> {
        std::vector<UNSIGNED_INTEGER_TYPE> outside_nodes;
        outside_nodes.reserve(std::floor(embedded_mesh->number_of_nodes() / 10.)); // Reserve 10% of the mesh size for outside nodes

        std::vector<BarycentricPoint> inside_nodes;
        inside_nodes.reserve(embedded_mesh->number_of_nodes());

        for (UNSIGNED_INTEGER_TYPE node_id = 0; node_id < embedded_mesh->number_of_nodes(); ++node_id) {
            BarycentricPoint b_point = barycentric_point(embedded_mesh->position(node_id));
            if (b_point.element_index < 0) {
                outside_nodes.emplace_back(node_id);
            } else {
                inside_nodes.emplace_back(b_point);
            }
        }

        // Add the mesh only if no nodes were found outside
        if (outside_nodes.empty()) {
            p_embedded_meshes[embedded_mesh] = inside_nodes;
        }

        return outside_nodes;
    }

    /**
     * Interpolate a field (scalar or vector field) from the container domain to the embedded mesh.
     *
     * This methods take an input field values (one value per node of the container mesh) and interpolate it onto
     * the values of the embedded nodes (outputs one value per node of the embedded mesh).
     *
     * @tparam Derived1 The matrix (Eigen) type of the input field values.
     * @tparam Derived2 The matrix (Eigen) type of the output field valuse.
     * @param embedded_mesh [INPUT] The embedded mesh on which the values will be interpolated. It must have been added
     *                              to the BarycentricContainer beforehand using the BarycentricContainer::add_embedded_mesh
     *                              method.
     * @param container_field_values [INPUT] The field values on every nodes of the container domain's mesh. This should be a
     *                                       matrix (Eigen) having one field value (scalar or vector) per rows. The number of
     *                                       rows should match the number of nodes of the container domain's mesh.
     * @param embedded_field_values [OUTPUT] The matrix (Eigen) where the interpolated field values should be written to.
     *                                       This should be a matrix (Eigen) having one field value (scalar or vector) per
     *                                       rows. The number of rows should match the number of nodes of the embedded mesh.
     */
    template <typename Derived1, typename Derived2>
    void interpolate_field(const BaseMesh & embedded_mesh,
                           const Eigen::MatrixBase<Derived1> & container_field_values,
                           Eigen::MatrixBase<Derived2> & embedded_field_values) const {

        // Make sure the embedded mesh has been registered.
        const auto mesh_entry = p_embedded_meshes.find(&embedded_mesh);
        if (mesh_entry == p_embedded_meshes.end()) {
            throw std::runtime_error("The embedded mesh hasn't been registered to the barycentric container. "
                                     "Make sure you called the add_embedded_mesh method before calling this one.");
        }

        // Make sure we have one field value per node of the container domain's mesh.
        if (static_cast<unsigned>(container_field_values.rows()) != p_container_domain->mesh().number_of_nodes()) {
            throw std::runtime_error("The number of rows of the input matrix must be the same as the number of nodes"
                                     "in the container domain's mesh.");
        }

        // Make sure we can write one field value per node of the embedded mesh.
        if (static_cast<unsigned>(embedded_field_values.rows()) != embedded_mesh.number_of_nodes()) {
            throw std::runtime_error("The number of rows of the output matrix must be the same as the number of nodes"
                                     "in the embedded mesh.");
        }

        if (container_field_values.cols() != embedded_field_values.cols()) {
            throw std::runtime_error("The number of columns of the input matrix (container field values) must be the "
                                     "same as the number of columns of the output matrix (embedded field values).");
        }

        // Get the barycentric coordinates of the embedded nodes
        const std::vector<BarycentricPoint> & barycentric_points = mesh_entry->second;
        if (embedded_mesh.number_of_nodes() != barycentric_points.size()) {
            throw std::runtime_error("The number of registered barycentric points does not match the number of nodes"
                                     "in the embedded mesh. Did the mesh size changed since it was added to the "
                                     "BarycentricContainer?");
        }

        // Interpolate field values
        for (UNSIGNED_INTEGER_TYPE node_id = 0; node_id < embedded_mesh.number_of_nodes(); ++node_id) {
            const auto & element_index = barycentric_points[node_id].element_index;
            const auto & local_coordinates = barycentric_points[node_id].local_coordinates;
            const auto & element_node_indices = p_container_domain->element_indices(element_index);

            const auto element = ContainerElement(p_container_domain->mesh().positions(element_node_indices));

            // Get the field values at the nodes of the containing element
            Eigen::Matrix<FLOATING_POINT_TYPE, ContainerElement::NumberOfNodesAtCompileTime, Derived1::ColsAtCompileTime>
                values (element.number_of_nodes(), container_field_values.cols());

            for (UNSIGNED_INTEGER_TYPE i = 0; i < element.number_of_nodes(); ++i) {
                values.row(i) = container_field_values.row(element_node_indices[i]);
            }

            // Interpolate the field value at the local position
            const auto interpolated_value = element.interpolate(local_coordinates, values);

            // Write back the interpolated value in the output embedded field
            embedded_field_values.row(node_id) = interpolated_value;
        }
    }

    /**
     * Interpolate a field (scalar or vector field) from the container domain to the embedded nodes.
     *
     * This methods take an input field values (one value per node of the container mesh) and interpolate it onto
     * the values of the embedded nodes (outputs one value per embedded node).
     *
     * @tparam Derived1 The matrix (Eigen) type of the embedded node positions.
     * @tparam Derived2 The matrix (Eigen) type of the input field values.
     * @tparam Derived2 The matrix (Eigen) type of the output field valuse.
     * @param embedded_mesh [INPUT] The embedded nodes on which the values will be interpolated.
     * @param container_field_values [INPUT] The field values on every nodes of the container domain's mesh. This should be a
     *                                       matrix (Eigen) having one field value (scalar or vector) per rows. The number of
     *                                       rows should match the number of nodes of the container domain's mesh.
     * @param embedded_field_values [OUTPUT] The matrix (Eigen) where the interpolated field values should be written to.
     *                                       This should be a matrix (Eigen) having one field value (scalar or vector) per
     *                                       rows. The number of rows should match the number of embedded nodes.
     */
    template <typename Derived1, typename Derived2, typename Derived3>
    void interpolate_field(const Eigen::MatrixBase<Derived1> & embedded_positions,
                           const Eigen::MatrixBase<Derived2> & container_field_values,
                           Eigen::MatrixBase<Derived3> & embedded_field_values) const {

        // Make sure we have one field value per node of the container domain's mesh.
        if (static_cast<unsigned>(container_field_values.rows()) != p_container_domain->mesh().number_of_nodes()) {
            throw std::runtime_error("The number of rows of the input matrix must be the same as the number of nodes"
                                     "in the container domain's mesh.");
        }

        // Make sure we can write one field value per node of the embedded mesh.
        if (static_cast<unsigned>(embedded_field_values.rows()) != embedded_positions.rows()) {
            throw std::runtime_error("The number of rows of the output matrix must be the same as the number of embedded nodes.");
        }

        if (container_field_values.cols() != embedded_field_values.cols()) {
            throw std::runtime_error("The number of columns of the input matrix (container field values) must be the "
                                     "same as the number of columns of the output matrix (embedded field values).");
        }

        // Get the barycentric coordinates of the embedded nodes
        const auto b_points = barycentric_points(embedded_positions);

        // Make sure every points are contained inside an element of the domain
        std::size_t number_of_outside_points = 0;
        for (const BarycentricPoint & b_point : b_points) {
            if (b_point.element_index < 0) {
                ++number_of_outside_points;
            }
        }

        if (number_of_outside_points > 0) {
            throw std::runtime_error("There is " + std::to_string(number_of_outside_points) +
                                     " points outside of the domain (ie no element is embedding them).");
        }

        // Interpolate field values
        const auto number_of_embedded_points = embedded_positions.rows();
        for (Eigen::Index node_id = 0; node_id < number_of_embedded_points; ++node_id) {
            const auto & element_index = b_points[node_id].element_index;
            const auto & local_coordinates = b_points[node_id].local_coordinates;
            const auto & element_node_indices = p_container_domain->element_indices(element_index);

            const auto element = ContainerElement(p_container_domain->mesh().positions(element_node_indices));

            // Get the field values at the nodes of the containing element
            Eigen::Matrix<typename Eigen::MatrixBase<Derived2>::Scalar,
                          ContainerElement::NumberOfNodesAtCompileTime,
                          Eigen::MatrixBase<Derived2>::ColsAtCompileTime>
                values (element.number_of_nodes(), container_field_values.cols());

            for (UNSIGNED_INTEGER_TYPE i = 0; i < element.number_of_nodes(); ++i) {
                values.row(i) = container_field_values.row(element_node_indices[i]);
            }

            // Interpolate the field value at the local position
            const auto interpolated_value = element.interpolate(local_coordinates, values);

            // Write back the interpolated value in the output embedded field
            embedded_field_values.row(node_id) = interpolated_value;
        }
    }
private:
    const Domain * p_container_domain;
    std::unordered_map<const BaseMesh *, std::vector<BarycentricPoint>> p_embedded_meshes;
    std::unique_ptr<HashGridT> p_hash_grid;
};

} // namespace caribou::topology