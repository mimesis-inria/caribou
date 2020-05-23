#pragma once

#include <Caribou/config.h>
#include <Caribou/Constants.h>
#include <Caribou/Geometry/Element.h>
#include <Eigen/Core>


namespace caribou::geometry {

/**
 * Base class for hexahedral elements.
 * @tparam Derived Derived hexahedral type (usually a leaf type, which is a type that is used by the user, eg Hexahedron).
 */
template<typename Derived>
struct BaseHexahedron : public Element<Derived> {
    // Types
    using Base = Element<Derived>;

    using LocalCoordinates = typename Base::LocalCoordinates;
    using WorldCoordinates = typename Base::WorldCoordinates;

    using GaussNode = typename Base::GaussNode;

    template <UNSIGNED_INTEGER_TYPE Dim>
    using Vector = typename Base::template Vector<Dim>;

    template <UNSIGNED_INTEGER_TYPE Rows, UNSIGNED_INTEGER_TYPE Cols>
    using Matrix = typename Base::template Matrix<Rows, Cols>;

    // Constants
    static constexpr auto CanonicalDimension = Base::CanonicalDimension;
    static constexpr auto Dimension = Base::Dimension;
    static constexpr auto NumberOfNodesAtCompileTime = Base::NumberOfNodesAtCompileTime;
    static constexpr auto NumberOfGaussNodesAtCompileTime = Base::NumberOfGaussNodesAtCompileTime;

    static_assert(Dimension == 3, "Hexahedrons can only be of dimension 3.");

    /** Default empty constructor */
    BaseHexahedron() = default;

    /** Constructor from an Eigen matrix containing the positions of the hexahedron's nodes */
    template<typename EigenType, REQUIRES(EigenType::RowsAtCompileTime == NumberOfNodesAtCompileTime)>
    explicit BaseHexahedron(Eigen::EigenBase<EigenType> & nodes) :p_nodes(nodes) {}

    /** Constructor from an Eigen matrix containing the positions of the hexahedron's nodes */
    template<typename EigenType, REQUIRES(EigenType::RowsAtCompileTime == NumberOfNodesAtCompileTime)>
    explicit BaseHexahedron(const Eigen::EigenBase<EigenType> & nodes) :p_nodes(nodes) {}

    /** Constructor from a serie of nodes. */
    template <
        typename ...Nodes,
        REQUIRES(NumberOfNodesAtCompileTime == sizeof...(Nodes)+1)
    >
    explicit BaseHexahedron(const WorldCoordinates & first_node, Nodes&&...remaining_nodes)
    {
        construct_from_nodes<0>(first_node, std::forward<Nodes>(remaining_nodes)...);
    }

    // Public methods common to all hexahedral types

    /**
     * Get the list of node indices of the faces.
     * \sa Element::boundary_elements_node_indices
     */
    inline auto faces() const {
        return self().get_boundary_elements_nodes();
    }

    /**
     * Extract the orthogonal frame of the element by computing the cross product of the unit vectors
     * from the center position to its projection on opposite faces.
     *
     * \sa BaseHexahedron::frame(const LocalCoordinates & local_point)
     *
     * \warning If the hexahedron isn't rectangular, the frame extracted by this function will be a rough
     *          approximation that could be far from the real solution, especially for strongly deformed or
     *          inverted hexahedrons.
     */
    inline auto frame() const -> Matrix<3, 3>
    {
        return frame(LocalCoordinates::Zeros());
    }

    /**
     * Extract the frame positioned at the given position (in local coordinates) on the hexa by computing the cross
     * product of the unit vectors from the given position its projection on opposite faces.
     *
     * This function will return a matrix of the form:
     * | ux vx wx |
     * | uy vy wy |
     * | uz vz wz |
     *
     * Where (ux, uy, uz), (vx, vy, vz) and (wx, wy, wz) are orthogonal unitary vectors representing
     * the u, v and w frames in the current hexa. If the hexa is rectangular and not rotated, this matrix is the
     * Identity matrix. If it is rectangular but rotated, rotating the hexa by the transposed of this frame should
     * align the u,v ,w axis to the x,y,z world frame (identity matrix).
     *
     * \warning If the hexahedron isn't rectangular, the frame extracted by this function will be a rough
     *          approximation that could be far from the real solution, especially for strongly deformed or
     *          inverted hexahedrons.
     */
    inline auto frame(const LocalCoordinates & local_point) const -> Matrix<3, 3>
    {
        // Position of the point inside the hexahedron where the frame should be computed
        const auto p = this->world_coordinates( local_point );

        // Project of the point on the quad facing the u axis
        const auto projected_on_u = this->world_coordinates({1, local_point[1],local_point[2]});

        // Project of the point on the quad facing the v axis
        const auto projected_on_v = this->world_coordinates({local_point[0], 1,local_point[2]});

        /* @todo(jnbrunet2000@gmail.com): select between the pairs of axis (center-to-u, center-to-v),
           (center-to-u, center-to-w) and (center-to-v, center-to-w) to find the best match (closer to orthogonal) */

        // Vector from the point to its projection on the quad facing the u axis
        const auto point_to_u = projected_on_u - p;

        // Vector from the point to its projection on the quad facing the v axis
        const auto point_to_v = projected_on_v - p;

        // The u-axis of the computed frame
        const auto u = point_to_u.normalized();

        // The v-axis of the computed frame
        auto v = point_to_v.normalized();

        // The w-axis of the computed frame
        const WorldCoordinates w = u.cross(v).normalized();

        // v-axis (recompute the v-axis in case u and v aren't orthogonal
        v = w.cross(u).normalized();

        Matrix<3, 3> m;
        m << u, v, w;

        return m;
    }

private:
    // Implementations
    friend struct Element<Derived>;
    [[nodiscard]]
    inline auto get_number_of_nodes() const {return NumberOfNodesAtCompileTime;}
    inline auto get_number_of_gauss_nodes() const {return NumberOfGaussNodesAtCompileTime;}
    inline auto get_node(const UNSIGNED_INTEGER_TYPE & index) const {return WorldCoordinates(p_nodes.row(index));};
    inline auto get_nodes() const -> const auto & {return p_nodes;};
    inline auto get_center() const {return Base::world_coordinates(LocalCoordinates(0, 0, 0));};
    inline auto get_number_of_boundary_elements() const -> UNSIGNED_INTEGER_TYPE {return 6;};
    inline auto get_contains_local(const LocalCoordinates & xi, const FLOATING_POINT_TYPE & eps) const -> bool {
        const auto & u = xi[0];
        const auto & v = xi[1];
        const auto & w = xi[2];
        return IN_CLOSED_INTERVAL(-1-eps, u, 1+eps) and
               IN_CLOSED_INTERVAL(-1-eps, v, 1+eps) and
               IN_CLOSED_INTERVAL(-1-eps, w, 1+eps);
    }

    auto self() -> Derived& { return *static_cast<Derived*>(this); }
    auto self() const -> const Derived& { return *static_cast<const Derived*>(this); }

    template <size_t index, typename ...Nodes, REQUIRES(sizeof...(Nodes) >= 1)>
    inline
    void construct_from_nodes(const WorldCoordinates & first_node, Nodes&&...remaining_nodes) {
        p_nodes.row(index) = first_node;
        construct_from_nodes<index+1>(std::forward<Nodes>(remaining_nodes)...);
    }

    template <size_t index>
    inline
    void construct_from_nodes(const WorldCoordinates & last_node) {
        p_nodes.row(index) = last_node;
    }
protected:
    Matrix<NumberOfNodesAtCompileTime, Dimension> p_nodes;
};

}