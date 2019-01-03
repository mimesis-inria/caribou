#ifndef CARIBOU_GEOMETRY_QUAD_H
#define CARIBOU_GEOMETRY_QUAD_H

#include <Caribou/config.h>
#include <Caribou/Algebra/Vector.h>
#include <Caribou/Geometry/Internal/BaseQuad.h>
#include <Caribou/Geometry/Node.h>
#include <Caribou/Geometry/Interpolation/Quad.h>

namespace caribou {
namespace geometry {

template <size_t Dim, typename Interpolation_ = interpolation::Quad4>
struct Quad : public internal::BaseQuad<Dim, Interpolation_, Quad<Dim, Interpolation_>>
{
};

template <size_t Dim>
struct Quad <Dim, interpolation::Quad4> : public internal::BaseQuad<Dim, interpolation::Quad4, Quad<Dim, interpolation::Quad4>>
{

    using NodeType = caribou::geometry::Node<Dim>;
    using Index = std::size_t ;

    constexpr Quad() {};

    constexpr Quad(const NodeType & p0, const NodeType & p1, const NodeType & p2, const NodeType & p3)
    : p_nodes {p0, p1, p2, p3}
    {

    }

    constexpr
    NodeType
    node(Index index) const
    {
        return p_nodes[index];
    }

    constexpr
    NodeType &
    node(Index index)
    {
        return p_nodes[index];
    }

private:
    std::array<NodeType, 4> p_nodes;
};

} // namespace geometry
} // namespace caribou
#endif //CARIBOU_GEOMETRY_QUAD_H
