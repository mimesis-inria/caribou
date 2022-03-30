#pragma once

#include <Caribou/config.h>
#include <Caribou/constants.h>
#include <Caribou/Geometry/Element.h>
#include <Eigen/Core>


namespace caribou::geometry {

template<typename Derived>
struct BaseTetrahedron : public Element<Derived> {
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

    static_assert(Dimension == 3, "Tetrahedrons can only be of dimension 3.");

    /** Default empty constructor */
    BaseTetrahedron() = default;

    /** Constructor from an Eigen matrix containing the positions of the tetra's nodes */
    template<typename EigenType, REQUIRES(EigenType::RowsAtCompileTime == NumberOfNodesAtCompileTime)>
    explicit BaseTetrahedron(Eigen::EigenBase<EigenType> & nodes) :p_nodes(nodes.derived().template cast<typename Base::Scalar>()) {}

    /** Constructor from an Eigen matrix containing the positions of the tetra's nodes */
    template<typename EigenType, REQUIRES(EigenType::RowsAtCompileTime == NumberOfNodesAtCompileTime)>
    explicit BaseTetrahedron(const Eigen::EigenBase<EigenType> & nodes) :p_nodes(nodes.derived().template cast<typename Base::Scalar>()) {}

    /** Constructor from a serie of nodes. */
    template <
        typename ...Nodes,
        REQUIRES(NumberOfNodesAtCompileTime == sizeof...(Nodes)+1)
    >
    explicit BaseTetrahedron(const WorldCoordinates & first_node, Nodes&&...remaining_nodes)
    {
        construct_from_nodes<0>(first_node, std::forward<Nodes>(remaining_nodes)...);
    }

    // Public methods common to all tetrahedral types

    /**
     * Get the list of node indices of the faces.
     * \sa Element::boundary_elements_node_indices
     */
    inline auto faces() const {
        return self().get_boundary_elements_nodes();
    }

    /**
     * Extract the orthogonal frame of the element by computing the cross product of the unit vectors
     * from the position of the first node to its projection on opposite faces.
     *
     * This function will return a matrix of the form:
     * | ux vx wx |
     * | uy vy wy |
     * | uz vz wz |
     *
     * Where (ux, uy, uz), (vx, vy, vz) and (wx, wy, wz) are orthogonal unitary vectors representing
     * the u, v and w frames in the current element. If the element is orthogonal and not rotated, this matrix is the
     * Identity matrix. If it is orthogonal but rotated, rotating the element by the transposed of this frame should
     * align the u,v ,w axis to the x,y,z world frame (identity matrix).
     *
     * \warning If the element isn't orthogonal, the frame extracted by this function will be a rough
     *          approximation that could be far from the real solution, especially for strongly deformed or
     *          inverted elements.
     */
    inline auto frame() const -> Matrix<3, 3>
    {
        // The u-axis of the computed frame
        const auto u = (this->node(1) - this->node(0)).normalized();

        // The v-axis of the computed frame
        auto v = (this->node(2) - this->node(0)).normalized();

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
    inline auto get_center() const {return Base::world_coordinates(LocalCoordinates({1/4., 1/4., 1/4.}));};
    inline auto get_number_of_boundary_elements() const -> UNSIGNED_INTEGER_TYPE {return 4;};
    inline auto get_contains_local(const LocalCoordinates & xi, const FLOATING_POINT_TYPE & eps) const -> bool {
        const auto & u = xi[0];
        const auto & v = xi[1];
        const auto & w = xi[2];
        return (u > -eps) and (v > -eps) and (w > -eps) and (1 - u - v - w > -eps);
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