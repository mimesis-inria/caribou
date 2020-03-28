#pragma once

#include <Caribou/config.h>
#include <Caribou/Constants.h>
#include <Caribou/Geometry/Element.h>
#include <Eigen/Core>

namespace caribou::geometry {

template<typename Derived>
struct BaseTriangle : public Element<Derived> {
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

    static_assert(Dimension == 2 or Dimension == 3, "Triangles can only be of dimension 2 or 3.");

    /** Default empty constructor */
    BaseTriangle() = default;

    /** Constructor from an array of floating point type */
    BaseTriangle(const FLOATING_POINT_TYPE(nodes)[NumberOfNodesAtCompileTime]) :p_nodes(nodes) {}

    /** Constructor from a serie of nodes. */
    template <
        typename ...Nodes,
        REQUIRES(NumberOfNodesAtCompileTime == sizeof...(Nodes)+1)
    >
    BaseTriangle(const WorldCoordinates & first_node, Nodes&&...remaining_nodes)
    {
        construct_from_nodes<0>(first_node, std::forward<Nodes>(remaining_nodes)...);
    }

private:
    // Implementations
    friend class Element<Derived>;
    [[nodiscard]]
    inline auto get_number_of_nodes() const {return NumberOfNodesAtCompileTime;}
    inline auto get_number_of_gauss_nodes() const {return NumberOfGaussNodesAtCompileTime;}
    inline auto get_node(const UNSIGNED_INTEGER_TYPE & index) const {return WorldCoordinates(p_nodes.row(index));};
    inline auto get_nodes() const -> const auto & {return p_nodes;};
    inline auto get_center() const {return Base::T(LocalCoordinates({1/3., 1/3.}));};
    inline auto get_number_of_boundary_elements() const -> UNSIGNED_INTEGER_TYPE {return 3;};

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

template<UNSIGNED_INTEGER_TYPE _Dimension, UNSIGNED_INTEGER_TYPE _Order = Linear>
struct Triangle;

}