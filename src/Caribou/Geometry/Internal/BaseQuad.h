#ifndef CARIBOU_GEOMETRY_INTERNAL_BASEQUAD_H
#define CARIBOU_GEOMETRY_INTERNAL_BASEQUAD_H

#include <Caribou/config.h>
#include <Caribou/Algebra/Vector.h>
#include <Caribou/Geometry/Interpolation/InterpolationElement.h>

namespace caribou {
namespace geometry {
namespace internal {

template<size_t Dim, typename CanonicalElementType, typename QuadType>
struct BaseQuad : caribou::geometry::interpolation::InterpolationElement<QuadType, CanonicalElementType>
{
    static_assert(Dim == 2 or Dim == 3, "Quad can only be made in dimension 2 or 3.");
};

template<typename CanonicalElementType, typename QuadType>
struct BaseQuad<2, CanonicalElementType, QuadType> : caribou::geometry::interpolation::InterpolationElement<QuadType, CanonicalElementType>
{
};

template<typename CanonicalElementType, typename QuadType>
struct BaseQuad<3, CanonicalElementType, QuadType> : caribou::geometry::interpolation::InterpolationElement<QuadType, CanonicalElementType>
{
};

} // namespace internal
} // namespace geometry
} // namespace caribou

#endif //CARIBOU_GEOMETRY_INTERNAL_BASEQUAD_H
