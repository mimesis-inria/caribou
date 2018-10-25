#ifndef CARIBOU_GEOMETRY_POINT_H
#define CARIBOU_GEOMETRY_POINT_H

#include <Caribou/Geometry/Entity.h>
#include <cassert>
#include <numeric>

namespace caribou
{
namespace geometry
{

/**
 * A point in space (independent of the space dimension).
 * @tparam Dim Dimension of the current space (default to 3D).
 * @tparam TReal Type of the floating values.
 */
template<size_t Dim, typename TReal>
class BasePoint : public Entity
{
public:
    typedef BasePoint<Dim, TReal> Self;
    typedef TReal Real;
    static constexpr size_t Dimension = Dim;

    BasePoint() = default;

    BasePoint(TReal const (&c)[Dim]) : Entity() {
        std::copy(std::begin(c), std::end(c), std::begin(coordinate));
    }

    BasePoint(const Self & p) : Entity() {
        std::copy(std::begin(p.coordinate), std::end(p.coordinate), std::begin(coordinate));
    }

    inline Self &operator=(const Self & p) {
        // check for self-assignment
        if(&p == this)
            return *this;

        std::copy(std::begin(p.coordinate), std::end(p.coordinate), std::begin(coordinate));

        return *this;
    }

    template <typename TOtherReal>
    inline Self &operator=(const std::initializer_list<TOtherReal> * l) {
        static_assert(l->size() == Dim, "Cannot initialized a Point of n dimension with a list of m dimension");

        std::copy(std::begin(*l), std::end(*l), std::begin(coordinate));

        return *this;
    }

    template<size_t OtherDim, typename TOtherReal>
    inline bool operator==(const BasePoint<OtherDim, TOtherReal> & p) const {
        return (
                std::equal(std::begin(coordinate), std::end(coordinate), std::begin(p.coordinate))
        );
    }

    template<size_t OtherDim, typename TOtherReal>
    inline bool operator!=(const BasePoint<OtherDim, TOtherReal> & p) const {
        return not (*this == p);
    }

    inline TReal & operator[] (std::size_t x) {
        return coordinate[x];
    }

    inline const TReal & operator[] (std::size_t x) const {
        return coordinate[x];
    }

    template<typename TOtherReal>
    inline TReal operator*(const BasePoint<Dim, TOtherReal> & other) const {
        return std::inner_product(std::begin(coordinate), std::end(coordinate), std::begin(other.coordinate), 0);
    }

    std::array<TReal, Dim> coordinate;
};

template<size_t Dim, typename TReal=float>
class Point : public BasePoint<Dim, TReal>
{
};

/**
 * A 1D point in space.
 */
template<typename TReal>
class Point<1, TReal>  : public BasePoint<1, TReal>
{
    using Self = Point<1, TReal>;
public:
    Point() = default;

    Point(TReal const (&coordinates)[1]) : BasePoint<1, TReal>(coordinates) {}

    explicit Point (const TReal & x) {
        Self::coordinate[0] = x;
    }

    inline const TReal & x () const { return Self::coordinate[0]; }

    template <typename TOtherReal>
    inline void set_x(const TOtherReal & x) {Self::coordinate[0] = static_cast<TReal> (x);}

    template<typename TOtherReal>
    inline TReal operator*(const Point<1, TOtherReal> & other) const  {
        return
            x() * (TReal)other.x()
        ;
    }

};

template<typename TReal=float>
using Point1D = Point<1, TReal>;

/**
 * A 2D point in space.
 */
template<typename TReal>
class Point<2, TReal>  : public BasePoint<2, TReal>
{
    using Self = Point<2, TReal>;
public:
    Point() = default;

    Point(TReal const (&coordinates)[2]) : BasePoint<2, TReal>(coordinates) {}

    Point (const TReal & x, const TReal & y) {
        Self::coordinate[0] = x;
        Self::coordinate[1] = y;
    }

    inline const TReal & x () const { return Self::coordinate[0]; }
    inline const TReal & y () const { return Self::coordinate[1]; }

    inline void set_x(const TReal & x) {Self::coordinate[0] = x;}
    inline void set_y(const TReal & y) {Self::coordinate[1] = y;}

    template<typename TOtherReal>
    inline TReal operator*(const Point<2, TOtherReal> & other) const  {
        return
            x() * (TReal)other.x()
            +
            y() * (TReal)other.y()
        ;
    }
};

template<typename TReal=float>
using Point2D = Point<2, TReal>;

/**
 * A 3D point in space.
 */
template<typename TReal>
class Point<3, TReal>  : public BasePoint<3, TReal>
{
    using Self = Point<3, TReal>;
public:
    Point() = default;

    Point(TReal const (&coordinates)[3]) : BasePoint<3, TReal>(coordinates) {}

    Point (const TReal & x, const TReal & y, const TReal & z) {
        Self::coordinate[0] = x;
        Self::coordinate[1] = y;
        Self::coordinate[2] = z;
    }

    inline const TReal & x () const { return Self::coordinate[0]; }
    inline const TReal & y () const { return Self::coordinate[1]; }
    inline const TReal & z () const { return Self::coordinate[2]; }

    inline void set_x(const TReal & x) {Self::coordinate[0] = x;}
    inline void set_y(const TReal & y) {Self::coordinate[1] = y;}
    inline void set_z(const TReal & z) {Self::coordinate[2] = z;}

    template<typename TOtherReal>
    inline TReal operator*(const Point<3, TOtherReal> & other) const {
        return
            x() * (TReal)other.x()
            +
            y() * (TReal)other.y()
            +
            z() * (TReal)other.z()
        ;
    }
};

template<typename TReal=float>
using Point3D = Point<3, TReal>;

template<size_t Dim, typename TReal=float>
Point<Dim, TReal> make_point(TReal const (&coordinates)[Dim]) {
    return Point<Dim, TReal>(coordinates);
}

template<typename... TReal>
auto make_point (TReal&&... coordinates) {
    return make_point<sizeof...(coordinates)>({ std::forward<TReal>(coordinates)... });
}


} // namespace geometry

} // namespace caribou
#endif //CARIBOU_GEOMETRY_POINT_H
