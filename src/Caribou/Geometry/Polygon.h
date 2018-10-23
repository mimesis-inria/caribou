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
template<int NumberOfSegments, typename TSegment, typename TData=BaseData>
class BasePolygon : public Entity<TData>
{
public:
    typedef BasePolygon<NumberOfSegments, TSegment, TData> Self;
    typedef TSegment SegmentType;
    typedef typename SegmentType::PointType PointType;
    static constexpr int Dimension = PointType::Dimension;

    static_assert(Dimension >= 2, "A polygon must be defined on space of dimension 2 or greater.");

    BasePolygon() = delete;

    BasePolygon(TSegment const (&segments)[NumberOfSegments], const TData data = TData()) : Entity<TData>(data) {
        assert(
                segments[0][0] == segments[NumberOfSegments-1][1] &&
                        "The segments of a polygon must form a closed loop (the first node of the first segment must be "
                                "equal to the last node of the last segment)."
        );

        for (int i = 0; i < NumberOfSegments; ++i) {
            nodes[i] = segments[i][0];
        }
    }

protected:
    std::array<TSegment, NumberOfSegments> nodes;
};


template<int NumberOfSegments, typename TSegment, typename TData=BaseData>
class Polygon : public BasePolygon<NumberOfSegments, TSegment, TData>
{};
} // namespace geometry

} // namespace caribou

#endif //CARIBOU_GEOMETRY_POLYGON_H
