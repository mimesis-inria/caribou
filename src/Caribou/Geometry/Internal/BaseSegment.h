#ifndef CARIBOU_GEOMETRY_INTERNAL_BASESEGMENT_H
#define CARIBOU_GEOMETRY_INTERNAL_BASESEGMENT_H

namespace caribou::geometry::internal {

template<size_t Dim, typename CanonicalElementType, typename SegmentType>
struct BaseSegment : public CanonicalElementType
{
    static_assert(Dim == 1 or Dim == 2 or Dim == 3, "Segment can only be made in dimension 1, 2 or 3.");
};

} // namespace caribou::geometry::internal

#endif //CARIBOU_GEOMETRY_INTERNAL_BASESEGMENT_H
