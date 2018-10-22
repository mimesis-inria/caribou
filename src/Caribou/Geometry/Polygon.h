#ifndef CARIBOU_GEOMETRY_POLYGON_H
#define CARIBOU_GEOMETRY_POLYGON_H


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
template<typename TSegment, typename TData=BaseData>
class Polygon : public Entity<TData>
{
public:
    Polygon() = delete;

protected:
    std::array<TPoint, 2> nodes;
};

} // namespace geometry

} // namespace caribou

#endif //CARIBOU_GEOMETRY_POLYGON_H
