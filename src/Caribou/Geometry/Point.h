#ifndef CARIBOU_GEOMETRY_POINT_H
#define CARIBOU_GEOMETRY_POINT_H

#include <Caribou/Geometry/Entity.h>
#include <cassert>

namespace caribou
{
namespace geometry
{

/**
 * A point in space (independent of the space dimension).
 * @tparam Dim Dimension of the current space (default to 3D).
 * @tparam TReal Type of the floating values.
 */
template<int Dim, typename TData, typename TReal>
class BasePoint : public Entity<TData>
{
public:
    typedef BasePoint<Dim, TData, TReal> Self;
    typedef TReal Real;
    typedef TData Data;

    BasePoint() = default;

    BasePoint(TReal const (&c)[Dim], const TData data = TData()) : Entity<TData>(data) {
        std::copy(std::begin(c), std::end(c), std::begin(coordinate));
    }

    BasePoint(const Self & p) : Entity<TData>(p.data) {
        std::copy(std::begin(p.coordinate), std::end(p.coordinate), std::begin(coordinate));
    }

    inline Self &operator=(const Self & p) {
        // check for self-assignment
        if(&p == this)
            return *this;

        std::copy(std::begin(p.coordinate), std::end(p.coordinate), std::begin(coordinate));
        this->data = p.data;

        return *this;
    }

    inline Self &operator=(const std::initializer_list<TReal> * l) {
        assert(l->size() == Dim && "Cannot initialized a Point of n dimension with a list of m dimension");
        std::copy(std::begin(*l), std::end(*l), std::begin(coordinate));

        return *this;
    }

    inline bool operator==(const Self & p) const {
        return (
                this->data == p.data &&
                std::equal(std::begin(coordinate), std::end(coordinate), std::begin(p.coordinate))
        );
    }

    inline bool operator!=(const Self & p) const {
        return not (*this == p);
    }

    inline TReal & operator[] (std::size_t x) {
        return coordinate[x];
    }

protected:
    std::array<TReal, Dim> coordinate;
};

///**
// * A 1D point in space.
// */
//template<typename Data=BaseData, typename Real=float>
//class Point1D : public Point<1, Data, Real>
//{
//    using Parent = Point<1, Data, Real>;
//public:
//    Point1D () = default;
//    explicit Point1D (const Real & x) {
//        Parent::coordinate[0] = x;
//    }
//
//    inline const Real & x () const { return Parent::coordinate[0]; }
//    inline void set_x(const Real & x) {Parent::coordinate[0] = x;}
//};
//
///**
// * A 2D point in space.
// */
//template<typename Data=BaseData, typename Real=float>
//class Point2D :  public Point<2, Data, Real>
//{
//    using Parent = Point<2, Data, Real>;
//public:
//    Point2D () = default;
//    Point2D (const Real & x, const Real & y) {
//        Parent::coordinate[0] = x;
//        Parent::coordinate[1] = y;
//    }
//
//    inline const Real & x () const { return Parent::coordinate[0]; }
//    inline const Real & y () const { return Parent::coordinate[1]; }
//
//    inline void set_x(const Real & x) {Parent::coordinate[0] = x;}
//    inline void set_y(const Real & y) {Parent::coordinate[1] = y;}
//};

template<int Dim, typename TData=BaseData, typename TReal=float>
class Point : public BasePoint<Dim, TData, TReal>
{
};

/**
 * A 1D point in space.
 */
template<typename TData, typename TReal>
class Point<1, TData, TReal>  : public BasePoint<1, TData, TReal>
{
    using Self = Point<1, TData, TReal>;
public:
    Point() = default;

    Point(TReal const (&coordinates)[1], const TData data = TData()) : BasePoint<1, TData, TReal>(coordinates, data) {}

    explicit Point (const TReal & x) {
        Self::coordinate[0] = x;
    }

    inline const TReal & x () const { return Self::coordinate[0]; }

    inline void set_x(const TReal & x) {Self::coordinate[0] = x;}
};

template<typename TData=BaseData, typename TReal=float>
using Point1D = Point<1, TData, TReal>;

/**
 * A 2D point in space.
 */
template<typename TData, typename TReal>
class Point<2, TData, TReal>  : public BasePoint<2, TData, TReal>
{
    using Self = Point<2, TData, TReal>;
public:
    Point() = default;

    Point(TReal const (&coordinates)[2], const TData data = TData()) : BasePoint<2, TData, TReal>(coordinates, data) {}

    Point (const TReal & x, const TReal & y) {
        Self::coordinate[0] = x;
        Self::coordinate[1] = y;
    }

    inline const TReal & x () const { return Self::coordinate[0]; }
    inline const TReal & y () const { return Self::coordinate[1]; }

    inline void set_x(const TReal & x) {Self::coordinate[0] = x;}
    inline void set_y(const TReal & y) {Self::coordinate[1] = y;}
};

template<typename TData=BaseData, typename TReal=float>
using Point2D = Point<2, TData, TReal>;

/**
 * A 3D point in space.
 */
template<typename TData, typename TReal>
class Point<3, TData, TReal>  : public BasePoint<3, TData, TReal>
{
    using Self = Point<3, TData, TReal>;
public:
    Point() = default;

    Point(TReal const (&coordinates)[3], const TData data = TData()) : BasePoint<3, TData, TReal>(coordinates, data) {}

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
};

template<typename TData=BaseData, typename TReal=float>
using Point3D = Point<3, TData, TReal>;

template<int Dim, typename TData=BaseData, typename TReal=float>
Point<Dim, TData, TReal> make_point(TReal const (&coordinates)[Dim], const TData & data = TData()) {
    return Point<Dim, TData, TReal>(coordinates, data);
}

template<typename... TReal>
auto make_point (TReal&&... coordinates) {
    return make_point<sizeof...(coordinates)>({ std::forward<TReal>(coordinates)... });
}


} // namespace geometry

} // namespace caribou
#endif //CARIBOU_GEOMETRY_POINT_H
