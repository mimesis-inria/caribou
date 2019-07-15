#ifndef CARIBOU_GEOMETRY_SEGMENT_H
#define CARIBOU_GEOMETRY_SEGMENT_H


#include <Caribou/config.h>
#include <Caribou/Algebra/Vector.h>
#include <Caribou/Geometry/Node.h>
#include <Caribou/Geometry/Interpolation/Segment.h>
#include <Caribou/Geometry/Internal/BaseSegment.h>

namespace caribou {
namespace geometry {

template <size_t Dim, typename CanonicalElementType = interpolation::Segment2>
struct Segment : public internal::BaseSegment<Dim, CanonicalElementType, Segment<Dim, CanonicalElementType>>
{
    static constexpr INTEGER_TYPE NumberOfNodes = CanonicalElementType::NumberOfNodes;
    using NodeType = caribou::geometry::Node<Dim>;
    using Index = std::size_t ;
    using Real = FLOATING_POINT_TYPE;

    using LocalCoordinates = algebra::Vector<1, Real>;
    using WorldCoordinates = algebra::Vector<Dim, Real>;

    static_assert(Dim == 1 or Dim == 2 or Dim == 3, "Only 1D, 2D and 3D segments are supported.");

    template <
        typename ...Nodes,
        REQUIRES(NumberOfNodes == sizeof...(Nodes)+1)
    >
    constexpr
    Segment(const NodeType & first_node, Nodes&&...remaining_nodes)
        : p_nodes {first_node, std::forward<Nodes>(remaining_nodes)...}
    {}

    inline constexpr bool
    operator==(const Segment<Dim, CanonicalElementType> & other) const noexcept
    {
        return std::equal(this->nodes().begin(), this->nodes().end(), other.nodes().begin(), [](const NodeType & n1, const NodeType & n2) -> bool {
            return n1 == n2;
        });
    }

    /** Get the Node at given index */
    constexpr
    const NodeType &
    node(Index index) const
    {
        return p_nodes[index];
    }

    /** Get the Node at given index */
    constexpr
    NodeType &
    node(Index index)
    {
        return p_nodes[index];
    }

    /** Get a reference to the set of nodes */
    inline
    const std::array<NodeType, NumberOfNodes> &
    nodes() const
    {
        return p_nodes;
    }

    /** Compute the center position **/
    auto
    center() const noexcept
    {
        return T({0});
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
    algebra::Matrix<Dim, 2, Real>
    jacobian (const LocalCoordinates & coordinates) const
    {
        return CanonicalElementType::Jacobian(coordinates, p_nodes);
    }

private:
    std::array<NodeType, NumberOfNodes> p_nodes;
};

} // namespace geometry
} // namespace caribou

#endif //CARIBOU_GEOMETRY_SEGMENT_H
