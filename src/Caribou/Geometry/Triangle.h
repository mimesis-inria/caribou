#ifndef CARIBOU_GEOMETRY_TRIANGLE_H
#define CARIBOU_GEOMETRY_TRIANGLE_H

#include <Caribou/config.h>
#include <Caribou/Algebra/Vector.h>
#include <Caribou/Geometry/Internal/BaseTriangle.h>
#include <Caribou/Geometry/Node.h>
#include <Caribou/Geometry/Interpolation/Triangle.h>

namespace caribou {
namespace geometry {

template <size_t Dim, typename Interpolation_ = interpolation::Triangle3>
struct Triangle : public internal::BaseTriangle<Dim, Interpolation_, Triangle<Dim, Interpolation_>>
{
};

template <size_t Dim>
struct Triangle <Dim, interpolation::Triangle3> : public internal::BaseTriangle<Dim, interpolation::Triangle3, Triangle<Dim, interpolation::Triangle3>>
{

    using NodeType = caribou::geometry::Node<Dim>;
    using Index = std::size_t ;

    constexpr Triangle() {};

    constexpr Triangle(const NodeType & p0, const NodeType & p1, const NodeType & p2)
            : p_nodes {p0, p1, p2}
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
    std::array<NodeType, 3> p_nodes;
};

} // namespace geometry
} // namespace caribou
#endif //CARIBOU_GEOMETRY_TRIANGLE_H
