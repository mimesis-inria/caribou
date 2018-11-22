#ifndef CARIBOU_GEOMETRY_POINT_H
#define CARIBOU_GEOMETRY_POINT_H

#include <cassert>
#include <numeric>
#include <cstddef>
#include <array>
#include <initializer_list>
#include <algorithm>

#include <Caribou/config.h>
#include <Caribou/Algebra/Vector.h>

namespace caribou
{
namespace geometry
{

/**
 * A point in space (independent of the space dimension).
 * @tparam Dim Dimension of the current space (default to 3D).
 * @tparam PointType The derived point class
 */
template<size_t Dim, typename PointType>
class BasePoint
{
public:
    static constexpr size_t Dimension = Dim;
    using VectorType = caribou::algebra::Vector<Dimension, FLOATING_POINT_TYPE>;
    using ValueType = typename VectorType::ValueType;

    BasePoint() = default;

    template <typename OtherValueType>
    BasePoint(std::initializer_list<OtherValueType> il) : coordinates(il) {}

    BasePoint(const PointType & other) : coordinates(other.coordinates) {}

    BasePoint(const VectorType & coordinates) : coordinates(coordinates) {}

    template <typename OtherValueType>
    BasePoint(OtherValueType const (&coordinates)[Dimension]) : coordinates(coordinates) {}

    inline PointType
    scale(VectorType s) const
    {
        VectorType scaled_coordinates;
        for (size_t i = 0; i < Dimension; ++i) {
            scaled_coordinates[i] = coordinates[i] * s[i];
        }

        return PointType(scaled_coordinates);

    }

    inline PointType
    &operator=(const PointType & other)
    {
        // check for self-assignment
        if(&other == this)
            return *this;

        coordinates = other.coordinates;

        return *this;
    }

    inline PointType
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
    operator==(const PointType & other) const
    { return coordinates == other.coordinates; }

    inline bool
    operator!=(const PointType & other) const
    { return !(*this == other); }


    VectorType coordinates;
};

template<size_t Dim>
class Point : public BasePoint<Dim, Point<Dim>>
{
};

/**
 * A 1D point in space.
 */
 template <>
class Point<1>  : public BasePoint<1, Point<1>>
{
    using Self = Point<1>;
public:
    static constexpr size_t Dimension = 1;
    using ValueType = typename BasePoint<Dimension, Self>::ValueType;
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
class Point<2>  : public BasePoint<2, Point<2>>
{
    using Self = Point<2>;
public:
    static constexpr size_t Dimension = 2;
    using ValueType = typename BasePoint<Dimension, Self>::ValueType;
    Point() = default;

    Point (const ValueType & x, const ValueType & y) {
        Self::coordinates[0] = x;
        Self::coordinates[1] = y;
    }

    Point(std::initializer_list<ValueType> il) : BasePoint<Dimension, Self>() {
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
class Point<3>  : public BasePoint<3, Point<3>>
{
    using Self = Point<3>;
public:
    static constexpr size_t Dimension = 3;
    using ValueType = typename BasePoint<Dimension, Self>::ValueType;
    Point() = default;

    Point (const ValueType & x, const ValueType & y, const ValueType & z) {
        Self::coordinates[0] = x;
        Self::coordinates[1] = y;
        Self::coordinates[2] = z;
    }

    using BasePoint<3, Point<3>>::BasePoint;

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
    return Point<sizeof...(args) +1> (coordinates);
}


} // namespace geometry

} // namespace caribou
#endif //CARIBOU_GEOMETRY_POINT_H
