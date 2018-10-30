#ifndef CARIBOU_GEOMETRY_POLYGON_H
#define CARIBOU_GEOMETRY_POLYGON_H
#include <cassert>

#include <Caribou/Geometry/Point.h>
#include <Caribou/Geometry/Segment.h>

namespace caribou
{
namespace geometry
{

/**
 * A polygon (of n nodes) in space (independent of the space dimension).
 */
template<size_t NNodes, size_t Dim>
class Polygon
{
public:
    static constexpr size_t Dimension = Dim;
    static constexpr size_t NumberOfNodes = NNodes;

    using PointType = Point<Dimension>;
    using VectorType = typename PointType::VectorType;
    using SegmentType = Segment<Dimension>;

    static_assert(Dimension >= 2, "A polygon must be defined on space of dimension 2 or greater.");
    static_assert(NNodes >= 3, "A polygon must have at least three nodes.");

    Polygon() = default;

    Polygon(PointType const (&n)[NumberOfNodes]) : nodes(n) {}
    Polygon(const std::array<PointType, NumberOfNodes> & n) : nodes(n) {}
    Polygon(const std::initializer_list<PointType > & il) {
        std::copy(std::begin(il), std::end(il), std::begin(nodes));
    }

    inline VectorType normal() const {
        auto v1 = segment(0).direction();
        auto v2 = segment(1).direction();

        return v1.cross(v2).unit();
    }

    inline SegmentType segment(size_t index) const
    {
        size_t i1 = index;
        size_t i2 = index+1;

        if (i1 >= Dimension)
            i1 = i1 % Dimension;

        if (i2 >= Dimension)
            i2 = i2 % Dimension;

        return make_segment(nodes[i1], nodes[i2]);
    }

    inline bool operator==(const Polygon<NumberOfNodes, Dimension> & other) const {
        return (
                std::equal(std::begin(nodes), std::end(nodes), std::begin(other.nodes), [](const PointType& a, const PointType &b) -> bool {
                    return a == b;
                })
        );
    }

    inline bool operator!=(const Polygon<NumberOfNodes, Dimension> & other) const {
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
template<size_t NSegments, size_t Dimension>
Polygon<NSegments, Dimension>
make_polygon_from_segments(Segment<Dimension> const (&segments)[NSegments])
{
    assert(
            segments[0][0] == segments[NSegments-1][1] &&
            "The segments of a polygon must form a closed loop (the first node of the first segment must be "
            "equal to the last node of the last segment)."
    );

    using PointType = Point<Dimension>;

    std::array<PointType, NSegments> points;
    for (size_t i = 0; i < NSegments; ++i)
        points[i] = segments[i][0];

    return Polygon<NSegments, Dimension>(points);
}

template<size_t Dimension, typename... Args>
auto make_polygon (const Segment<Dimension> & arg, Args&&... args) {
    return make_polygon_from_segments<sizeof...(args)+1>({ arg, std::forward<Args>(args)... });
}

//
// Tool functions to create a polygon with a list of position vectors
//
template<size_t Dimension, typename... Args>
auto make_polygon (const Point<Dimension> & arg, Args&&... args)
{
    return Polygon<sizeof...(args)+1, Dimension>({ arg, std::forward<Args>(args)... });
}

template<size_t Dimension, size_t NumberOfNodes, typename ValueType>
auto make_polygon(ValueType const (&arg)[NumberOfNodes][Dimension])
{
    Polygon<NumberOfNodes, Dimension> p;
    for (size_t i = 0; i < NumberOfNodes; ++i) {
        p[i] = arg[i];
    }
    return p;
}

} // namespace geometry

} // namespace caribou

#endif //CARIBOU_GEOMETRY_POLYGON_H
