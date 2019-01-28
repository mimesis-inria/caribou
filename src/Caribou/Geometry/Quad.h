#ifndef CARIBOU_GEOMETRY_QUAD_H
#define CARIBOU_GEOMETRY_QUAD_H

#include <Caribou/config.h>
#include <Caribou/Algebra/Vector.h>
#include <Caribou/Geometry/Node.h>
#include <Caribou/Geometry/Interpolation/Quad.h>

namespace caribou {
namespace geometry {

template <size_t Dim, typename CanonicalElementType = interpolation::Quad4>
struct Quad : public CanonicalElementType
{
    static constexpr INTEGER_TYPE NumberOfNodes = CanonicalElementType::NumberOfNodes;
    using NodeType = caribou::geometry::Node<Dim>;
    using Index = std::size_t ;
    using Real = FLOATING_POINT_TYPE;

    static_assert(Dim == 2 or Dim == 3, "Only 2D and 3D quads are supported.");

    template <
            typename ...Nodes,
            REQUIRES(NumberOfNodes == sizeof...(Nodes)),
            REQUIRES(std::conjunction_v<std::is_same<NodeType, Nodes>...>)
    >
    constexpr
    Quad(Nodes&&...remaining_nodes)
    : p_nodes {std::forward<Nodes>(remaining_nodes)...}
    {}

    constexpr
    Quad(const std::array<NodeType, NumberOfNodes> & nodes)
            : p_nodes(nodes)
    {}

    const NodeType &
    node(Index index) const
    {
        return p_nodes[index];
    }

    NodeType &
    node(Index index)
    {
        return p_nodes[index];
    }

    /**
     * Compute the transformation of a local position {u,v} to its world position {x,y,z}
     */
    inline
    NodeType
    T(const Real & u, const Real & v) const
    {
        return CanonicalElementType::interpolate_at_local_position({u,v}, p_nodes);
    }

    /** Compute the jacobian matrix evaluated at local position {u,v} */
    algebra::Matrix<Dim, 2, Real>
    jacobian (const Real & u, const Real & v) const
    {
        return CanonicalElementType::Jacobian({u,v}, p_nodes);
    }

private:
    std::array<NodeType, NumberOfNodes> p_nodes;
};

} // namespace geometry
} // namespace caribou
#endif //CARIBOU_GEOMETRY_QUAD_H
