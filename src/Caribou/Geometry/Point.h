#ifndef CARIBOU_GEOMETRY_POINT_H
#define CARIBOU_GEOMETRY_POINT_H

#include <Caribou/Geometry/Entity.h>
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
 * @tparam TVector Type of the vector position. Default to caribou::algebra::Vector<Dim>>
 */
template<size_t Dim, typename TVector = caribou::algebra::Vector<Dim>>
class BasePoint : public Entity
{
public:
    using VectorType = TVector;
    static constexpr size_t Dimension = Dim;

    BasePoint() = default;

    template <typename ValueType>
    BasePoint(std::initializer_list<ValueType> il) : Entity() , coordinates(il) {}

    template<typename OtherVectorType>
    BasePoint(const BasePoint<Dim, OtherVectorType> & other) : Entity(), coordinates(other.coordinates) {}

    template<typename OtherVectorType>
    BasePoint(const OtherVectorType & coordinates) : Entity(), coordinates(coordinates) {}

    template<typename OtherVectorType>
    inline BasePoint<Dim, TVector>
    &operator=(const BasePoint<Dim, OtherVectorType> & other)
    {
        // check for self-assignment
        if(&other == this)
            return *this;

        coordinates = other.coordinates;

        return *this;
    }

    inline BasePoint<Dim, TVector>
    &operator=(const TVector & coordinates)
    {
        this->coordinates = coordinates;

        return *this;
    }

    template <typename ValueType>
    inline ValueType &
    operator[] (std::size_t x)
    {
        return coordinates[x];
    }

    template <typename ValueType>
    inline const ValueType &
    operator[] (std::size_t x) const
    {
        return coordinates[x];
    }

    // TODO(jnbrunet2000@gmail.com): This should not be mandatory for the vector template. We should add it using the detection idiom
    using ValueType = typename VectorType::value_type;
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


    template<typename OtherPoint>
    inline bool
    operator==(const OtherPoint & other) const
    { return coordinates == other.coordinates; }

    template<typename OtherPoint>
    inline bool
    operator!=(const OtherPoint & other) const
    { return !(*this == other); }


    VectorType coordinates;
};

template<size_t Dim, typename TVector = caribou::algebra::Vector<Dim>>
class Point : public BasePoint<Dim, TVector>
{
};

/**
 * A 1D point in space.
 */
template<typename TVector>
class Point<1, TVector>  : public BasePoint<1, TVector>
{
    using Self = Point<1, TVector>;
    // TODO(jnbrunet2000@gmail.com): This should not be mandatory for the vector template. We should add it using the detection idiom
    using ValueType = typename BasePoint<1, TVector>::ValueType;
public:
    Point() = default;

    template <typename ValueType>
    explicit Point (const ValueType & x) {
        Self::coordinates[0] = x;
    }

    inline const ValueType & x () const { return Self::coordinates[0]; }
    inline ValueType & x () { return Self::coordinates[0]; }
};

/**
 * A 2D point in space.
 */
template<typename TVector>
class Point<2, TVector>  : public BasePoint<2, TVector>
{
    using Self = Point<2, TVector>;
    // TODO(jnbrunet2000@gmail.com): This should not be mandatory for the vector template. We should add it using the detection idiom
    using ValueType = typename BasePoint<2, TVector>::ValueType;
public:
    Point() = default;

    template <typename ValueType>
    Point (const ValueType & x, const ValueType & y) {
        Self::coordinates[0] = x;
        Self::coordinates[1] = y;
    }

    template <typename ValueType>
    Point(std::initializer_list<ValueType> il) : BasePoint<2, TVector>() {
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
template<typename TVector>
class Point<3, TVector>  : public BasePoint<3, TVector>
{
    using Self = Point<3, TVector>;
    // TODO(jnbrunet2000@gmail.com): This should not be mandatory for the vector template. We should add it using the detection idiom
    using ValueType = typename BasePoint<3, TVector>::ValueType;
public:
    Point() = default;

    template <typename ValueType>
    Point (const ValueType & x, const ValueType & y, const ValueType & z) {
        Self::coordinates[0] = x;
        Self::coordinates[1] = y;
        Self::coordinates[2] = z;
    }

    template <typename ValueType>
    Point(std::initializer_list<ValueType> il) : BasePoint<3, TVector>() {
        std::copy(std::begin(il), std::end(il), std::begin(Self::coordinates));
    }

    template <typename ValueType>
    Point(ValueType const (&coordinates)[3]) : BasePoint<3, TVector>() {
        Self::coordinates = coordinates;
    }

    inline const ValueType & x () const { return Self::coordinates[0]; }
    inline ValueType & x () { return Self::coordinates[0]; }

    inline const ValueType & y () const { return Self::coordinates[1]; }
    inline ValueType & y () { return Self::coordinates[1]; }

    inline const ValueType & z () const { return Self::coordinates[2]; }
    inline ValueType & z () { return Self::coordinates[2]; }
};

template<typename TVector = caribou::algebra::Vector<1>>
using Point1D = Point<1, TVector>;

template<typename TVector = caribou::algebra::Vector<2>>
using Point2D = Point<2, TVector>;

template<typename TVector = caribou::algebra::Vector<3>>
using Point3D = Point<3, TVector>;

template<typename ValueType, typename... Args>
auto make_point (ValueType value, Args&&... args) {
    caribou::algebra::Vector<sizeof...(args) + 1, ValueType> coordinates ({ value, std::forward<Args>(args)... });
    return coordinates;
}


} // namespace geometry

} // namespace caribou
#endif //CARIBOU_GEOMETRY_POINT_H
