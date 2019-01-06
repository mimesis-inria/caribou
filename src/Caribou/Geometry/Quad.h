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
    static constexpr INTEGER_TYPE NumberOfNodes = Interpolation_::NumberOfNodes;
    using NodeType = caribou::geometry::Node<Dim>;
    using Index = std::size_t ;

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
    const NodeType &
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
    std::array<NodeType, Interpolation_::NumberOfNodes> p_nodes;
};

} // namespace geometry
} // namespace caribou
#endif //CARIBOU_GEOMETRY_QUAD_H
