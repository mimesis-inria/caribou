#ifndef CARIBOU_GEOMETRY_INTERNAL_BASEQUAD_H
#define CARIBOU_GEOMETRY_INTERNAL_BASEQUAD_H

#include <Caribou/config.h>
#include <Caribou/Algebra/Vector.h>

namespace caribou {
namespace geometry {
namespace internal {

template<size_t Dim, typename Interpolation, typename QuadType>
struct BaseQuad : Interpolation
{
    static_assert(Dim == 2 or Dim == 3, "Quad can only be made in dimension 2 or 3.");
};

template<typename Interpolation, typename QuadType>
struct BaseQuad<2, Interpolation, QuadType> : Interpolation
{
};

template<typename Interpolation, typename QuadType>
struct BaseQuad<3, Interpolation, QuadType> : Interpolation
{
};

} // namespace internal
} // namespace geometry
} // namespace caribou

#endif //CARIBOU_GEOMETRY_INTERNAL_BASEQUAD_H
