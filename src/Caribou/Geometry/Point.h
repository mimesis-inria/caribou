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
template<int Dim=3, typename TData=BaseData, typename TReal=float>
class Point : public Entity<TData>
{
public:
    typedef TReal Real;
    typedef TData Data;

    Point() = default;
    Point(std::initializer_list<Real> l) : coordinate(l) {}
    Point(const Point<Dim, Data, Real> & p) {
        std::copy(std::begin(p.coordinate), std::end(p.coordinate), std::begin(coordinate));
        this->data = p.data;
    }

    inline Point<Dim, Data, Real> &operator=(const Point<Dim, Data, Real> & p) {
        // check for self-assignment
        if(&p == this)
            return *this;

        std::copy(std::begin(p.coordinate), std::end(p.coordinate), std::begin(coordinate));
        this->data = p.data;

        return *this;
    }

    inline Point<Dim, Data, Real> &operator=(const std::initializer_list<Real> * l) {
        assert(l->size() == Dim && "Cannot initialized a Point of n dimension with a list of m dimension");
        std::copy(std::begin(*l), std::end(*l), std::begin(coordinate));

        return *this;
    }

    inline bool operator==(const Point<Dim, Data, Real> & p) const {
        return (
                this->data == p.data &&
                std::equal(std::begin(coordinate), std::end(coordinate), std::begin(p.coordinate))
        );
    }

    inline bool operator!=(const Point<Dim, Data, Real> & p) const {
        return not (*this == p);
    }

    inline Real & operator[] (std::size_t x) {
        return coordinate[x];
    }

protected:
    std::array<Real, Dim> coordinate;
};

/**
 * A 1D point in space.
 */
template<typename Data=BaseData, typename Real=float>
class Point1D : public Point<1, Data, Real>
{
    using Parent = Point<1, Data, Real>;
public:
    Point1D () = default;
    explicit Point1D (const Real & x) {
        Parent::coordinate[0] = x;
    }

    inline const Real & x () const { return Parent::coordinate[0]; }
    inline void set_x(const Real & x) {Parent::coordinate[0] = x;}
};

/**
 * A 2D point in space.
 */
template<typename Data=BaseData, typename Real=float>
class Point2D :  public Point<2, Data, Real>
{
    using Parent = Point<2, Data, Real>;
public:
    Point2D () = default;
    Point2D (const Real & x, const Real & y) {
        Parent::coordinate[0] = x;
        Parent::coordinate[1] = y;
    }

    inline const Real & x () const { return Parent::coordinate[0]; }
    inline const Real & y () const { return Parent::coordinate[1]; }

    inline void set_x(const Real & x) {Parent::coordinate[0] = x;}
    inline void set_y(const Real & y) {Parent::coordinate[1] = y;}
};

/**
 * A 3D point in space.
 */
template<typename Data=BaseData, typename Real=float>
class Point3D :  public Point<3, Data, Real>
{
    using Parent = Point<3, Data, Real>;
public:
    Point3D () = default;
    Point3D (const Real & x, const Real & y, const Real & z) {
        Parent::coordinate[0] = x;
        Parent::coordinate[1] = y;
        Parent::coordinate[2] = z;
    }

    inline const Real & x () const { return Parent::coordinate[0]; }
    inline const Real & y () const { return Parent::coordinate[1]; }
    inline const Real & z () const { return Parent::coordinate[2]; }

    inline void set_x(const Real & x) {Parent::coordinate[0] = x;}
    inline void set_y(const Real & y) {Parent::coordinate[1] = y;}
    inline void set_z(const Real & z) {Parent::coordinate[2] = z;}
};

} // namespace geometry

} // namespace caribou
#endif //CARIBOU_GEOMETRY_POINT_H
