#ifndef CARIBOU_GEOMETRY_HEXAHEDRON_H
#define CARIBOU_GEOMETRY_HEXAHEDRON_H

#include <Caribou/config.h>
#include <Caribou/Algebra/Vector.h>
#include <Caribou/Geometry/Node.h>
#include <Caribou/Geometry/Interpolation/Hexahedron.h>

namespace caribou {
namespace geometry {

template <typename CanonicalElementType = interpolation::Hexahedron8>
struct Hexahedron : public CanonicalElementType
{
    static constexpr INTEGER_TYPE NumberOfNodes = CanonicalElementType::NumberOfNodes;
    using NodeType = caribou::geometry::Node<3>;
    using Index = std::size_t ;
    using Real = FLOATING_POINT_TYPE;

    template <
            typename ...Nodes,
            REQUIRES(NumberOfNodes == sizeof...(Nodes)),
            REQUIRES(std::conjunction_v<std::is_same<NodeType, Nodes>...>)
    >
    constexpr
    Hexahedron(Nodes&&...remaining_nodes)
    : p_nodes {std::forward<Nodes>(remaining_nodes)...}
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

    /** Compute the jacobian matrix evaluated at local position {u,v} */
    algebra::Matrix<3, 3, Real>
    jacobian (const Real & u, const Real & v, const Real & w) const
    {
        return CanonicalElementType::Jacobian({u,v, w}, p_nodes);
    }

private:
    std::array<NodeType, NumberOfNodes> p_nodes;
};

} // namespace geometry
} // namespace caribou
#endif //CARIBOU_GEOMETRY_HEXAHEDRON_H
