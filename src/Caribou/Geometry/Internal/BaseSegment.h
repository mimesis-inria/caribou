#ifndef CARIBOU_GEOMETRY_INTERNAL_BASESEGMENT_H
#define CARIBOU_GEOMETRY_INTERNAL_BASESEGMENT_H

#include <cmath>

#include <Caribou/config.h>
#include <Caribou/Algebra/Vector.h>

namespace caribou {
namespace geometry {
namespace internal {

template<size_t Dim, typename CanonicalElementType, typename SegmentType>
struct BaseSegment : public CanonicalElementType
{
    static_assert(Dim == 1 or Dim == 2 or Dim == 3, "Segment can only be made in dimension 1, 2 or 3.");
};

} // namespace internal
} // namespace geometry
} // namespace caribou

#endif //CARIBOU_GEOMETRY_INTERNAL_BASESEGMENT_H
