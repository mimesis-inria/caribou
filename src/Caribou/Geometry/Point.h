#ifndef CARIBOU_GEOMETRY_POINT_H
#define CARIBOU_GEOMETRY_POINT_H

#include <Caribou/Algebra/Vector.h>
#include <cassert>
#include <numeric>
#include <cstddef>
#include <array>
#include <initializer_list>
#include <algorithm>

namespace caribou
{
namespace geometry
{

/**
 * A point in space (independent of the space dimension).
 * @tparam Dim Dimension of the current space (default to 3D).
 */
template<size_t Dim>
class BasePoint
{
public:
    using VectorType = caribou::algebra::Vector<Dim>;
    using ValueType = typename VectorType::ValueType;
    static constexpr size_t Dimension = Dim;

    BasePoint() = default;

    template <typename ValueType>
    BasePoint(std::initializer_list<ValueType> il) : coordinates(il) {}

    BasePoint(const BasePoint<Dim> & other) : coordinates(other.coordinates) {}

    BasePoint(const VectorType & coordinates) : coordinates(coordinates) {}

    inline BasePoint<Dimension>
    scale(VectorType s) const
    {
        VectorType scaled_coordinates;
        for (size_t i = 0; i < Dimension; ++i) {
            scaled_coordinates[i] = coordinates[i] * s[i];
        }

        return BasePoint<Dimension>(scaled_coordinates);

    }

    inline BasePoint<Dim>
    &operator=(const BasePoint<Dim> & other)
    {
        // check for self-assignment
        if(&other == this)
            return *this;

        coordinates = other.coordinates;

        return *this;
    }

    inline BasePoint<Dim>
    &operator=(const VectorType & coordinates)
    {
        this->coordinates = coordinates;

        return *this;
    }

    inline ValueType &
    operator[] (std::size_t x)
    {
        return coordinates[x];
    }

    inline const ValueType &
    operator[] (std::size_t x) const
    {
        return coordinates[x];
    }

    inline bool
    operator==(const BasePoint<Dimension> & other) const
    { return coordinates == other.coordinates; }

    inline bool
    operator!=(const BasePoint<Dimension> & other) const
    { return !(*this == other); }


    VectorType coordinates;
};

template<size_t Dim>
class Point : public BasePoint<Dim>
{
};

/**
 * A 1D point in space.
 */
 template <>
class Point<1>  : public BasePoint<1>
{
    using Self = Point<1>;
public:
    using ValueType = typename BasePoint<1>::ValueType;
    Point() = default;

    explicit Point (const ValueType & x) {
        Self::coordinates[0] = x;
    }

    inline const ValueType & x () const { return Self::coordinates[0]; }
    inline ValueType & x () { return Self::coordinates[0]; }
};

/**
 * A 2D point in space.
 */
template<>
class Point<2>  : public BasePoint<2>
{
    using Self = Point<2>;
public:
    using ValueType = typename BasePoint<2>::ValueType;
    Point() = default;

    Point (const ValueType & x, const ValueType & y) {
        Self::coordinates[0] = x;
        Self::coordinates[1] = y;
    }

    Point(std::initializer_list<ValueType> il) : BasePoint<2>() {
        std::copy(std::begin(il), std::end(il), std::begin(Self::coordinates));
    }

    inline const ValueType & x () const { return Self::coordinates[0]; }
    inline ValueType & x () { return Self::coordinates[0]; }

    inline const ValueType & y () const { return Self::coordinates[1]; }
    inline ValueType & y () { return Self::coordinates[1]; }
};

/**
 * A 3D point in space.
 */
template<>
class Point<3>  : public BasePoint<3>
{
    using Self = Point<3>;
public:
    using ValueType = typename BasePoint<3>::ValueType;
    Point() = default;

    Point (const ValueType & x, const ValueType & y, const ValueType & z) {
        Self::coordinates[0] = x;
        Self::coordinates[1] = y;
        Self::coordinates[2] = z;
    }

    template <typename OtherValueType>
    Point(std::initializer_list<OtherValueType> il) : BasePoint<3>() {
        std::copy(std::begin(il), std::end(il), std::begin(Self::coordinates));
    }

    template <typename OtherValueType>
    Point(OtherValueType const (&coordinates)[3]) : BasePoint<3>() {
        Self::coordinates = coordinates;
    }

    inline const ValueType & x () const { return Self::coordinates[0]; }
    inline ValueType & x () { return Self::coordinates[0]; }

    inline const ValueType & y () const { return Self::coordinates[1]; }
    inline ValueType & y () { return Self::coordinates[1]; }

    inline const ValueType & z () const { return Self::coordinates[2]; }
    inline ValueType & z () { return Self::coordinates[2]; }
};

using Point1D = Point<1>;
using Point2D = Point<2>;
using Point3D = Point<3>;

template<typename ValueType, typename... Args>
auto make_point (ValueType value, Args&&... args) {
    caribou::algebra::Vector<sizeof...(args) + 1, ValueType> coordinates ({ value, std::forward<Args>(args)... });
    return BasePoint<sizeof...(args) +1> (coordinates);
}


} // namespace geometry

} // namespace caribou
#endif //CARIBOU_GEOMETRY_POINT_H
