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
 * A polygon (of n segments) in space (independent of the space dimension).
 */
template<size_t NSegments, typename TSegment, typename TData=BaseData>
class Polygon : public Entity<TData>
{
public:
    typedef TSegment SegmentType;
    typedef typename SegmentType::PointType PointType;
    static constexpr size_t Dimension = PointType::Dimension;
    static constexpr size_t NumberOfSegments = NSegments;
    static constexpr size_t NumberOfNodes = NSegments;

    static_assert(Dimension >= 2, "A polygon must be defined on space of dimension 2 or greater.");
    static_assert(NSegments >= 3, "A polygon must have at least three segments.");

    Polygon() = default;

    Polygon(TSegment const (&s)[NumberOfSegments], const TData data = TData()) : Entity<TData>(data) {
        assert(
                s[0][0] == s[NumberOfSegments-1][1] &&
                        "The segments of a polygon must form a closed loop (the first node of the first segment must be "
                                "equal to the last node of the last segment)."
        );

        for (size_t i = 0; i < NumberOfSegments; ++i) {
            nodes[i] = s[i][0];
        }
    }

    Polygon(PointType const (&n)[NumberOfNodes], const TData data = TData()) : Entity<TData>(data) {
        for (size_t i = 0; i < NumberOfNodes; ++i)
            nodes[i] = n[i];
    }

    template<typename TOtherSegment, typename TOtherData>
    inline bool operator==(const Polygon<NSegments, TOtherSegment, TOtherData> & other) const {
        using OtherPointType = typename TOtherSegment::PointType;
        return (
                this->data == other.data &&
                std::equal(std::begin(nodes), std::end(nodes), std::begin(other.nodes), [](const PointType& a, const OtherPointType &b) -> bool {
                    return a == b;
                })
        );
    }

    template<typename TOtherSegment, typename TOtherData>
    inline bool operator!=(const Polygon<NSegments, TOtherSegment, TOtherData> & other) const {
        return not (*this == other);
    }

    inline PointType & operator[] (std::size_t x) {
        return nodes[x];
    }

    inline const PointType & operator[] (std::size_t x) const {
        return nodes[x];
    }

protected:
    std::array<PointType, NumberOfSegments> nodes;
};

//
// Tool functions to create a polygon with a list of segments
//
template<size_t NSegments, typename TSegment, typename TData=BaseData>
Polygon<NSegments, TSegment, TData>
make_polygon_from_segments(TSegment const (&segments)[NSegments], const TData &data = TData())
{
    return Polygon<NSegments, TSegment, TData>(segments, data);
}

template<class TPoint, typename TData, typename... Args>
auto make_polygon (Segment<TPoint, TData> & arg, Args&&... args) {
    return make_polygon_from_segments<sizeof...(args)+1>({ arg, std::forward<Args>(args)... });
}


//
// Tool functions to create a polygon with a list of points
//
template<size_t NPoints,
        typename TPoint,
        typename TData=BaseData
>
Polygon<NPoints, Segment<TPoint, BaseData>, TData>
make_polygon_from_points(TPoint const (&points)[NPoints], const TData &data = TData())
{
    return Polygon<NPoints, Segment<TPoint, BaseData>, TData>(points, data);
}

template<size_t Dim, typename TData, typename TReal, typename... Args>
auto make_polygon (Point<Dim, TData, TReal> & arg, Args&&... args) {
    return make_polygon_from_points<sizeof...(args)+1>({ arg, std::forward<Args>(args)... });
}

} // namespace geometry

} // namespace caribou

#endif //CARIBOU_GEOMETRY_POLYGON_H
