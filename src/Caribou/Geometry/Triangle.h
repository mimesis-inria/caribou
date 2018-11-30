#ifndef CARIBOU_GEOMETRY_TRIANGLE_H
#define CARIBOU_GEOMETRY_TRIANGLE_H

#include <Caribou/Geometry/Polygon.h>
#include <Caribou/Geometry/Segment.h>

namespace caribou
{
namespace geometry
{

/** A triangle is a polygon of three nodes. **/
template<size_t Dimension>
class Triangle : public Polygon<3, Dimension>
{
public:
    using Base = Polygon<3, Dimension>;
    using PointType = Point<Dimension>;
    using VectorType = typename PointType::VectorType;
    using SegmentType = Segment<Dimension>;
    using Float = FLOATING_POINT_TYPE;

    /** Copy constructor from another 3 points polygon. */
    Triangle(const Polygon<3, Dimension> & other) {
        std::copy(std::begin(other.nodes), std::end(other.nodes), std::begin(this->nodes));
    }

    /** Compute the surface area **/
    inline Float area() const noexcept {
        const VectorType v1 = Base::segment(0).direction();
        const VectorType v2 = Base::segment(1).direction();

        return v1.cross(v2).length() / 2.;
    }

    /** Compute the surface area **/
    inline PointType center() const noexcept {
        const PointType & p1 = Base::node(0);
        const PointType & p2 = Base::node(1);
        const PointType & p3 = Base::node(2);

        return (p1 + p2 + p3)/ (Float) 3.0;
    }

};

/**
 * Create a triangle from three points.
 *
 * Example:
 * \code{.cpp}
 * auto p1 = make_point(0,0,0);
 * auto p2 = make_point(1,1,1);
 * auto p3 = make_point(0,1,0);
 *
 * auto t = make_triangle(p1, p2, p3);
 * \endcode
 * @return
 */
template<size_t Dimension>
Triangle<Dimension>
make_triangle(const Point<Dimension> & p1, const Point<Dimension> & p2, const Point<Dimension> & p3) {
    return make_polygon(p1, p2, p3);
}

/** Create a triangle from three points of type PointType. */
template<size_t Dimension, typename PointType>
Triangle<Dimension>
make_triangle(const PointType & p1, const PointType& p2, const PointType & p3) {
    const Point<Dimension> point1 = p1;
    const Point<Dimension> point2 = p2;
    const Point<Dimension> point3 = p3;

    return make_triangle(point1, point2, point3);
}

/**
 * Create a triangle from a list of lists.
 *
 * Example:
 * \code{.cpp}
 * auto t = make_triangle(
 *   {{0, 0, 0}, {1, 1, 1}, {0, 1, 0}}
 * );
 * \endcode
 * @return
 */
template<size_t Dimension, typename ValueType>
auto
make_triangle(ValueType const (&arg)[3][Dimension])
{
    return make_polygon(arg);
}

/**
 * Create a triangle from three lists.
 *
 * Example:
 * \code{.cpp}
 * auto t = make_triangle(
 *   {0, 0, 0}, {1, 1, 1}, {0, 1, 0}
 * );
 * \endcode
 * @return
 */
template<size_t Dimension, typename ValueType>
auto
make_triangle(const ValueType (&arg1)[Dimension], const ValueType (&arg2)[Dimension], const ValueType (&arg3)[Dimension])
{
    ValueType nodes[3][Dimension];

    for (size_t i = 0; i < Dimension; ++i) {
        nodes[0][i] = arg1[i];
        nodes[1][i] = arg2[i];
        nodes[2][i] = arg3[i];
    }

    return make_polygon(nodes);
}

} // namespace geometry

} // namespace caribou

#endif //CARIBOU_GEOMETRY_TRIANGLE_H
