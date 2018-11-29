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
class BasePoint : public caribou::algebra::Vector<Dim, FLOATING_POINT_TYPE>
{
public:
    static constexpr size_t Dimension = Dim;
    using VectorType = caribou::algebra::Vector<Dimension, FLOATING_POINT_TYPE>;
    using ValueType = typename VectorType::ValueType;

    using VectorType::VectorType;

    BasePoint() : VectorType() {}
    BasePoint(bool initialize_to_zero) : VectorType(initialize_to_zero) {}
    BasePoint(const VectorType & v) : VectorType(v) {}
    BasePoint(const PointType & p) {
        for (size_t i = 0; i < Dimension; ++i)
            (*this)[i] = static_cast<ValueType> (p[i]);
    }

    /**
     * Return a copy of this point translated by the vector t
     */
    inline PointType
    translated(const VectorType & t) const
    {
        return PointType(*this).translate(t);
    }


    /**
     * Translate this point by the vector t
     */
    inline PointType &
    translate(const VectorType & t)
    {
        for (size_t i = 0; i < Dimension; ++i) {
            (*this)[i] += t[i];
        }

        return static_cast<PointType&>(*this);
    }

    /**
     * Return a copy of this point where each component are scaled by the corresponding scaling factor in s
     */
    inline PointType
    scaled(const VectorType & s) const
    {
        PointType scaled_coordinates;
        for (size_t i = 0; i < Dimension; ++i) {
            scaled_coordinates[i] *=  (*this)[i] * s[i];
        }

        return scaled_coordinates;
    }

    /**
     * Return a copy of this point with its coordinates scaled by the scalar s.
     */
    inline PointType
    scaled(const ValueType & s) const
    {
        return (*this) * s;
    }

    /**
     * Scale each component of this point's coordinates by the corresponding scaling factor in s
     */
    inline PointType &
    scale(const VectorType & s)
    {
        for (size_t i = 0; i < Dimension; ++i) {
            (*this)[i] *= s[i];
        }

        return static_cast<PointType&>(*this);
    }

    /**
     * Scale this point's coordinate by the scalar s.
     */
    inline PointType &
    scale(const ValueType & s)
    {
        for (size_t i = 0; i < Dimension; ++i) {
            (*this)[i] *= s;
        }
        return static_cast<PointType&>(*this);
    }
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
    using Base = BasePoint<1, Point<1>>;

public:
    static constexpr size_t Dimension = 1;
    using ValueType = typename Base::ValueType;

    using Base::Base; // Import constructors

    /** Assignment operator **/
    template<typename OtherValueType>
    Point &
    operator = (const caribou::algebra::Vector<Dimension, OtherValueType> & v)
    {
        for (size_t i = 0; i < Dimension; ++i)
            (*this)[i] = v[i];
        return (*this);
    }

    inline const ValueType & x () const { return this->at(0); }
    inline ValueType & x () { return this->at(0); }
};

/**
 * A 2D point in space.
 */
template<>
class Point<2>  : public BasePoint<2, Point<2>>
{
    using Self = Point<2>;
    using Base = BasePoint<2, Point<2>>;

public:
    static constexpr size_t Dimension = 2;
    using ValueType = typename BasePoint<Dimension, Self>::ValueType;

    using Base::Base; // Import constructors

    /** Assignment operator **/
    template<typename OtherValueType>
    Point &
    operator = (const caribou::algebra::Vector<Dimension, OtherValueType> & v)
    {
        for (size_t i = 0; i < Dimension; ++i)
            (*this)[i] = v[i];
        return (*this);
    }

    inline const ValueType & x () const { return this->at(0); }
    inline ValueType & x () { return this->at(0); }

    inline const ValueType & y () const { return this->at(1); }
    inline ValueType & y () { return this->at(1); }
};

/**
 * A 3D point in space.
 */
template<>
class Point<3>  : public BasePoint<3, Point<3>>
{
    using Self = Point<3>;
    using Base = BasePoint<3, Point<3>>;

public:
    static constexpr size_t Dimension = 3;
    using ValueType = typename BasePoint<Dimension, Self>::ValueType;

    using Base::Base;

    /** Assignment operator **/
    template<typename OtherValueType>
    Point &
    operator = (const caribou::algebra::Vector<Dimension, OtherValueType> & v)
    {
        for (size_t i = 0; i < Dimension; ++i)
            (*this)[i] = v[i];
        return (*this);
    }

    inline const ValueType & x () const { return this->at(0); }
    inline ValueType & x () { return this->at(0); }

    inline const ValueType & y () const { return this->at(1); }
    inline ValueType & y () { return this->at(1); }

    inline const ValueType & z () const { return this->at(2); }
    inline ValueType & z () { return this->at(2); }
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
