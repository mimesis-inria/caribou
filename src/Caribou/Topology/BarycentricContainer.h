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
 * The BarycentricContainer class allows to embed nodes into a container domain.
 * The embedded nodes do not need to match the nodes of the container domain. However, an embedded node must
 * be found in exactly one of the element of the container, or at the boundary between elements (for example, an
 * embedded node lying on a face between two elements of the container is valid). If an embedded node
 * is lying outside of the container domain (ie is not located inside any of the container
 * elements), the BarycentricContainer will ignore it. The list of ignored nodes can be retrieved.
 *
 * @tparam Domain The type of the container domain.
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
        BarycentricPoint() : element_index(-1), local_coordinates(LocalCoordinates::Zero()) {}
        BarycentricPoint(const ElementIndex & index, const LocalCoordinates & coordinates)
        : element_index(index), local_coordinates(coordinates) {}

        ElementIndex element_index;
        LocalCoordinates local_coordinates;
    };

    /** Default constructor is not permitted */
    BarycentricContainer() = delete;

    /**
     * Construct the container from the given domain. This will create an HashGrid class instance to be able
     * to quickly retrieve the container element of a given world position. The size of the HashGrid cells will
     * be set to the mean size of the container elements.
     *
     * This constructs the BarycentricContainer with a set of embedded points where the barycentric coordinates
     * will be computed.
     *
     * @tparam Derived NXD Eigen matrix representing the D dimensional coordinates of the N embedded positions.
     * @param container_domain The mesh domain that will contain the embedded meshes.
     * @param embedded_points The positions (in world coordinates) embedded in the container mesh for which the barycentric
     *                        points have to be found.
     */
    template <typename Derived>
    BarycentricContainer(const Domain * container_domain, const Eigen::MatrixBase<Derived> & embedded_points)
    : BarycentricContainer(container_domain) {
        // Set the embedded points
        set_embedded_points(embedded_points);
    }

    /**
     * Construct the container from the given domain. This will create an HashGrid class instance to be able
     * to quickly retrieve the container element of a given world position. The size of the HashGrid cells will
     * be set to the mean size of the container elements.
     *
     * @param container_domain The mesh domain that will contain the embedded nodes.
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
        // List of candidate elements that could contain the point
        const auto candidate_element_indices = p_hash_grid->get(p);


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
     * Get the list of closest elements to a point and its barycentric coordinates within these elements.
     */
    auto closest_elements(const WorldCoordinates & p) const -> std::vector<BarycentricPoint> {
        // List of candidate elements that could contain the point
        const auto candidate_element_indices = p_hash_grid->get(p);

        std::vector<BarycentricPoint> closest_elements;
        closest_elements.reserve(candidate_element_indices.size());
        for (const auto & element_index : candidate_element_indices) {
            const ContainerElement e = p_container_domain->element(element_index);
            const LocalCoordinates local_coordinates = e.local_coordinates(p);
            closest_elements.emplace_back(static_cast<ElementIndex>(element_index), local_coordinates);
        }

        return closest_elements;
    }

    /**
     * Interpolate a field (scalar or vector field) from the container domain to the embedded nodes.
     *
     * This methods take an input field values (one value per node of the container mesh) and interpolate it onto
     * the values of the embedded nodes (outputs one value per embedded node). Embedded nodes that are found outside of the
     * containing domain (ie, no elements containing the node can be found) are ignored and their interpolated values won't
     * be computed.
     *
     * @tparam Derived1 The matrix (Eigen) type of the input field values.
     * @tparam Derived2 The matrix (Eigen) type of the output field valuse.
     * @param container_field_values [INPUT] The field values on every nodes of the container domain. This should be a
     *                                       matrix (Eigen) having one field value (scalar or vector) per rows. The
     *                                       number of rows should match the number of nodes of the container domain.
     * @param embedded_field_values [OUTPUT] The matrix (Eigen) where the interpolated field values should be written to.
     *                                       This should be a matrix (Eigen) having one field value (scalar or vector)
     *                                       per rows. The number of rows should match the number of embedded nodes.
     */
    template <typename Derived1, typename Derived2>
    void interpolate(const Eigen::MatrixBase<Derived1> & container_field_values,
                           Eigen::MatrixBase<Derived2> & embedded_field_values) const {

        // Make sure we have one field value per node of the container domain's mesh.
        if (static_cast<unsigned>(container_field_values.rows()) != p_container_domain->mesh().number_of_nodes()) {
            std::stringstream ss;
            ss << "The number of rows of the input matrix must be the same as the number of nodes in the container domain. ";
            ss << "The number of rows is " << container_field_values.rows() << " and the number of nodes is ";
            ss << p_container_domain->mesh().number_of_nodes();
            throw std::runtime_error(ss.str());
        }

        // Make sure we can write one field value per node of the embedded mesh.
        if (static_cast<unsigned>(embedded_field_values.rows()) != p_barycentric_points.size()) {
            throw std::runtime_error("The number of rows of the output matrix must be the same as the number of embedded nodes.");
        }

        if (container_field_values.cols() != embedded_field_values.cols()) {
            throw std::runtime_error("The number of columns of the input matrix (container field values) must be the "
                                     "same as the number of columns of the output matrix (embedded field values).");
        }

        // Interpolate field values
        const auto number_of_embedded_points = static_cast<signed >(p_barycentric_points.size());
        for (Eigen::Index node_id = 0; node_id < number_of_embedded_points; ++node_id) {
            const auto & element_index = p_barycentric_points[node_id].element_index;
            const auto & local_coordinates = p_barycentric_points[node_id].local_coordinates;

            if (element_index < 0) {
                continue;
            }

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
            if constexpr (Eigen::MatrixBase<Derived2>::ColsAtCompileTime == 1) {
                embedded_field_values[node_id] = interpolated_value;
            } else {
                embedded_field_values.row(node_id) = interpolated_value;
            }
        }
    }

    /**
     * Barycentric points of the embedded nodes. See BarycentricContainer::BarycentricPoint for more details.
     */
    auto barycentric_points () const -> const std::vector<BarycentricPoint> & {
        return p_barycentric_points;
    }

    /**
     * Indices of nodes found outside of the domain (ie, nodes that were not found inside any container elements of the
     * domain).
     */
    [[nodiscard]]
    auto outside_nodes () const -> const std::vector<UNSIGNED_INTEGER_TYPE> & {
        return p_outside_nodes;
    }

private:
    /**
     * Set the barycentric points from a set of positions embedded inside the container domain.
     * @tparam Derived NXD Eigen matrix representing the D dimensional coordinates of the N embedded positions.
     * @param embedded_points The positions (in world coordinates) embedded in the container mesh for which the barycentric
     *                        points have to be found.
     * \note When an embedded position lie completely outside the container mesh (i.e. it is not contained inside any
     *       container element), its index will be added to the list of outside nodes and it will be ignored by future
     *       interpolation calls.
     */
    template <typename Derived>
    void set_embedded_points(const Eigen::MatrixBase<Derived> & embedded_points) {
        static_assert(Eigen::MatrixBase<Derived>::ColsAtCompileTime == Dimension or Eigen::MatrixBase<Derived>::ColsAtCompileTime == Eigen::Dynamic);

        if (Eigen::MatrixBase<Derived>::ColsAtCompileTime == Eigen::Dynamic and embedded_points.cols() != Dimension) {
            throw std::runtime_error("Trying to get the barycentric coordinates of " +
                                     std::to_string(embedded_points.cols()) + "D points from a " +
                                     std::to_string(Dimension) + "D domain");
        }

        p_outside_nodes.resize(0);
        p_outside_nodes.reserve(std::floor(embedded_points.rows() / 10.)); // Reserve 10% of the mesh size for outside nodes

        p_barycentric_points.resize(0);
        p_barycentric_points.reserve(embedded_points.rows());

        const auto number_of_embedded_points = embedded_points.rows();
        for (Eigen::Index node_id = 0; node_id < number_of_embedded_points; ++node_id) {
            const BarycentricPoint bp = barycentric_point(embedded_points.row(node_id).template cast<typename WorldCoordinates::Scalar>());
            if (bp.element_index < 0) {
                p_outside_nodes.emplace_back(node_id);
            }

            p_barycentric_points.emplace_back(bp);
        }
    }


    const Domain * p_container_domain;
    std::vector<BarycentricPoint> p_barycentric_points;
    std::vector<UNSIGNED_INTEGER_TYPE> p_outside_nodes;
    std::unique_ptr<HashGridT> p_hash_grid;
};

} // namespace caribou::topology