#ifndef CARIBOU_GEOMETRY_INTERNAL_BASETRIANGLE_H
#define CARIBOU_GEOMETRY_INTERNAL_BASETRIANGLE_H

#include <cmath>

#include <Caribou/config.h>
#include <Caribou/Algebra/Vector.h>

namespace caribou {
namespace geometry {
namespace internal {

template<size_t Dim, typename CanonicalElementType, typename TriangleType>
struct BaseTriangle : public CanonicalElementType
{
    static_assert(Dim == 2 or Dim == 3, "Triangle can only be made in dimension 2 or 3.");
};

template<typename CanonicalElementType, typename TriangleType>
struct BaseTriangle<2, CanonicalElementType, TriangleType> : public CanonicalElementType
{

private:
    const TriangleType &self () const
    {
        return static_cast<const TriangleType &>(*this);
    }
};

template<typename CanonicalElementType, typename TriangleType>
struct BaseTriangle<3, CanonicalElementType, TriangleType> : public CanonicalElementType
{
    inline auto normal() const noexcept {
        const auto v1 = self().node(1) - self().node(0);
        const auto v2 = self().node(2) - self().node(1);

        return v1.cross(v2).unit();
    }

private:
    const TriangleType &self () const
    {
        return static_cast<const TriangleType &>(*this);
    }
};

} // namespace internal
} // namespace geometry
} // namespace caribou

#endif //CARIBOU_GEOMETRY_INTERNAL_BASETRIANGLE_H
