#ifndef CARIBOU_GEOMETRY_POLYGON_H
#define CARIBOU_GEOMETRY_POLYGON_H
#include <cassert>

#include <Caribou/Geometry/Entity.h>
#include <Caribou/Geometry/Point.h>
#include <Caribou/Geometry/Segment.h>

namespace caribou
{
namespace geometry
{

/**
 * A polygon (of n nodes) in space (independent of the space dimension).
 */
template<size_t NNodes, typename TPoint>
class Polygon : public Entity
{
public:
    typedef TPoint PointType;
    static constexpr size_t Dimension = PointType::Dimension;
    static constexpr size_t NumberOfSegments = NNodes;
    static constexpr size_t NumberOfNodes = NNodes;

    static_assert(Dimension >= 2, "A polygon must be defined on space of dimension 2 or greater.");
    static_assert(NNodes >= 3, "A polygon must have at least three nodes.");

    Polygon() = default;

    Polygon(PointType const (&n)[NumberOfNodes]) : Entity(), nodes(n) {}
    Polygon(const std::array<PointType, NumberOfNodes> & n) : Entity(), nodes(n) {}
    Polygon(const std::initializer_list<PointType > & il) : Entity() {
        std::copy(std::begin(il), std::end(il), std::begin(nodes));
    }

    template<typename TOtherPoint>
    inline bool operator==(const Polygon<NNodes, TOtherPoint> & other) const {
        return (
                std::equal(std::begin(nodes), std::end(nodes), std::begin(other.nodes), [](const PointType& a, const TOtherPoint &b) -> bool {
                    return a == b;
                })
        );
    }

    template<typename TOtherPoint>
    inline bool operator!=(const Polygon<NNodes, TOtherPoint> & other) const {
        return not (*this == other);
    }

    inline PointType & operator[] (std::size_t x) {
        return nodes[x];
    }

    inline const PointType & operator[] (std::size_t x) const {
        return nodes[x];
    }

    std::array<PointType, NumberOfNodes> nodes;
};

//
// Tool functions to create a polygon with a list of segments
//
template<size_t NSegments, typename TPoint>
Polygon<NSegments, TPoint>
make_polygon_from_segments(Segment<TPoint> const (&segments)[NSegments])
{
    assert(
            segments[0][0] == segments[NSegments-1][1] &&
            "The segments of a polygon must form a closed loop (the first node of the first segment must be "
            "equal to the last node of the last segment)."
    );

    std::array<TPoint, NSegments> points;
    for (size_t i = 0; i < NSegments; ++i)
        points[i] = segments[i][0];

    return Polygon<NSegments, TPoint>(points);
}

template<class TPoint, typename... Args>
auto make_polygon (Segment<TPoint> & arg, Args&&... args) {
    return make_polygon_from_segments<sizeof...(args)+1>({ arg, std::forward<Args>(args)... });
}

//
// Tool functions to create a polygon with a list of points
//
template<size_t Dim, typename TReal, typename... Args>
auto make_polygon (Point<Dim, TReal> & arg, Args&&... args)
{
    return Polygon<sizeof...(args)+1, Point<Dim, TReal>>({ arg, std::forward<Args>(args)... });
}

template<size_t Dimension, size_t NumberOfNodes, typename ValueType>
auto make_polygon(ValueType const (&arg)[NumberOfNodes][Dimension])
{
    Polygon<NumberOfNodes, Point<Dimension, caribou::algebra::Vector<Dimension,ValueType>>> p;
    for (size_t i = 0; i < NumberOfNodes; ++i) {
        p[i] = arg[i];
    }
    return p;
}

} // namespace geometry

} // namespace caribou

#endif //CARIBOU_GEOMETRY_POLYGON_H
