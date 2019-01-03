#ifndef CARIBOU_GEOMETRY_INTERNAL_BASETRIANGLE_H
#define CARIBOU_GEOMETRY_INTERNAL_BASETRIANGLE_H

#include <cmath>

#include <Caribou/config.h>
#include <Caribou/Algebra/Vector.h>

namespace caribou {
namespace geometry {
namespace internal {

template<size_t Dim, typename Interpolation, typename TriangleType>
struct BaseTriangle : Interpolation
{
    static_assert(Dim == 2 or Dim == 3, "Triangle can only be made in dimension 2 or 3.");
};

template<typename Interpolation, typename TriangleType>
struct BaseTriangle<2, Interpolation, TriangleType> : Interpolation
{
    /** Compute the surface area **/
    inline auto area() const noexcept {
        const auto n1 = self().node(0);
        const auto n2 = self().node(1);
        const auto n3 = self().node(2);

        using ValueType = decltype(n1[0]);

        const caribou::algebra::Matrix M ({
                {n1[0], n2[0], n3[0]},
                {n1[1], n2[1], n3[1]},
                {(ValueType) 1., (ValueType) 1., (ValueType) 1.}
        });

        return 1/2. * std::abs(M.determinant());
    }

    /** Compute the surface area **/
    inline auto center() const noexcept {
        const auto & p1 = self().node(0);
        const auto & p2 = self().node(1);
        const auto & p3 = self().node(2);

        return (p1 + p2 + p3)/ 3.0;
    }

private:
    const TriangleType &self () const
    {
        return static_cast<const TriangleType &>(*this);
    }
};

template<typename Interpolation, typename TriangleType>
struct BaseTriangle<3, Interpolation, TriangleType> : Interpolation
{
    /** Compute the surface area **/
    inline auto area() const noexcept {
        const auto v1 = self().node(1) - self().node(0);
        const auto v2 = self().node(2) - self().node(1);

        return v1.cross(v2).length() / 2.;
    }

    /** Compute the surface area **/
    inline auto center() const noexcept {
        const auto & p1 = self().node(0);
        const auto & p2 = self().node(1);
        const auto & p3 = self().node(2);

        return (p1 + p2 + p3)/ 3.0;
    }

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
