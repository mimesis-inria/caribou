#ifndef CARIBOU_GEOMETRY_INTERNAL_BASESEGMENT_H
#define CARIBOU_GEOMETRY_INTERNAL_BASESEGMENT_H

#include <Eigen/Core>

namespace caribou::geometry::internal {

template<size_t Dim, typename CanonicalElementType, typename SegmentType>
struct BaseSegment : public CanonicalElementType
{
    static_assert(Dim == 1 or Dim == 2 or Dim == 3, "Segment can only be made in dimension 1, 2 or 3.");

    template<int nNodes>
    using NodesContainer = Eigen::Matrix<FLOATING_POINT_TYPE, nNodes, Dim, Eigen::RowMajor>;
};


template<typename CanonicalElementType, typename SegmentType>
struct BaseSegment<1, CanonicalElementType, SegmentType> : public CanonicalElementType
{
    template<int nNodes>
    using NodesContainer = Eigen::Matrix<FLOATING_POINT_TYPE, nNodes, 1, Eigen::ColMajor>;
};
} // namespace caribou::geometry::internal

#endif //CARIBOU_GEOMETRY_INTERNAL_BASESEGMENT_H
