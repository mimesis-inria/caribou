#ifndef CARIBOU_GEOMETRY_QUAD_H
#define CARIBOU_GEOMETRY_QUAD_H

#include <Eigen/Core>

#include <Caribou/config.h>
#include <Caribou/Traits.h>
#include <Caribou/Geometry/Interpolation/Quad.h>

namespace caribou::geometry {

template <size_t Dim, typename CanonicalElementType = interpolation::Quad4>
struct Quad : public CanonicalElementType
{
    static constexpr INTEGER_TYPE NumberOfNodes = CanonicalElementType::NumberOfNodes;

    using LocalCoordinates = typename CanonicalElementType::LocalCoordinates;
    using WorldCoordinates = Eigen::Matrix<FLOATING_POINT_TYPE, Dim, 1>;

    template<int nRows, int nColumns, int Options=Eigen::RowMajor>
    using Matrix = Eigen::Matrix<FLOATING_POINT_TYPE, nRows, nColumns, Options>;

    static_assert(Dim == 2 or Dim == 3, "Only 2D and 3D quads are supported.");

    template <
        typename ...Nodes,
        REQUIRES(NumberOfNodes == sizeof...(Nodes)+1)
    >
    Quad(const WorldCoordinates & first_node, Nodes&&...remaining_nodes)
    {
        construct_from_nodes<0>(first_node, std::forward<Nodes>(remaining_nodes)...);
    }

    Quad(const Matrix<NumberOfNodes, Dim> & nodes)
    : p_nodes(nodes)
    {}

    inline
    auto
    node(UNSIGNED_INTEGER_TYPE index) const
    {
        return p_nodes.row(index).transpose();
    }

    inline
    auto
    node(UNSIGNED_INTEGER_TYPE index)
    {
        return p_nodes.row(index).transpose();
    }

    /** Get a reference to the set of nodes */
    inline
    const auto &
    nodes() const
    {
        return p_nodes;
    }

    /** Compute the center position **/
    auto
    center() const noexcept
    {
        return T(LocalCoordinates({0., 0.}));
    }

    /**
     * Compute the transformation of a local position {u} to its world position {x,y,z}
     */
    inline
    WorldCoordinates
    T(const LocalCoordinates & coordinates) const
    {
        return CanonicalElementType::interpolate_at_local_position(coordinates, nodes());
    }

    /** Compute the jacobian matrix evaluated at local position {u,v} */
    Matrix<Dim, 2>
    jacobian (const LocalCoordinates & coordinates) const
    {
        return CanonicalElementType::Jacobian(coordinates, p_nodes);
    }

private:
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


private:
    Matrix<NumberOfNodes, Dim> p_nodes;
};

} // namespace caribou::geometry
#endif //CARIBOU_GEOMETRY_QUAD_H
