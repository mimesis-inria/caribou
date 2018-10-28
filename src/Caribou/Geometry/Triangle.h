#ifndef CARIBOU_GEOMETRY_TRIANGLE_H
#define CARIBOU_GEOMETRY_TRIANGLE_H

#include <Caribou/Geometry/Polygon.h>

namespace caribou
{
namespace geometry
{

/** A triangle is an alias to a polygon of three nodes. **/
template<typename VectorType>
using Triangle = Polygon<3, VectorType>;

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
template<typename TPoint>
Triangle<typename TPoint::VectorType> make_triangle(const TPoint & p1, const TPoint & p2, const TPoint & p3) {
    return make_polygon(p1, p2, p3);
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
auto make_triangle(ValueType const (&arg)[3][Dimension])
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
auto make_triangle(const ValueType (&arg1)[Dimension], const ValueType (&arg2)[Dimension], const ValueType (&arg3)[Dimension])
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
