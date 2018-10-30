#ifndef CARIBOU_GEOMETRY_QUAD_H
#define CARIBOU_GEOMETRY_QUAD_H

#include <Caribou/Geometry/Polygon.h>

namespace caribou
{
namespace geometry
{

/** A quad is an alias to a polygon of four nodes. **/
template<typename VectorType>
using Quad = Polygon<4, VectorType>;

/**
 * Create a quad from four points.
 *
 * Example:
 * \code{.cpp}
 * auto p1 = make_point(1,-1,0);
 * auto p2 = make_point(-1,-1,0);
 * auto p3 = make_point(1,1,0);
 * auto p4 = make_point(-1,1,0);
 *
 * auto q = make_quad(p1, p2, p3, p4);
 * \endcode
 * @return
 */
template<typename TPoint>
Quad<typename TPoint::VectorType>
make_quad(const TPoint & p1, const TPoint & p2, const TPoint & p3, const TPoint & p4) {
    return make_polygon(p1, p2, p3, p4);
}

/**
 * Create a quad from a list of lists.
 *
 * Example:
 * \code{.cpp}
 * auto q = make_quad(
 *   {{1, -1, 0}, {-1, -1, 0}, {1, 1, 0}, {-1, 1, 0}}
 * );
 * \endcode
 * @return
 */
template<size_t Dimension, typename ValueType>
auto
make_quad(ValueType const (&arg)[4][Dimension])
{
    return make_polygon(arg);
}

/**
 * Create a quad from four lists.
 *
 * Example:
 * \code{.cpp}
 * auto q = make_quad(
 *   {1, -1, 0}, {-1, -1, 0}, {1, 1, 0}, {-1, 1, 0}
 * );
 * \endcode
 * @return
 */
template<size_t Dimension, typename ValueType>
auto
make_quad(const ValueType (&arg1)[Dimension], const ValueType (&arg2)[Dimension], const ValueType (&arg3)[Dimension], const ValueType (&arg4)[Dimension])
{
    ValueType nodes[4][Dimension];

    for (size_t i = 0; i < Dimension; ++i) {
        nodes[0][i] = arg1[i];
        nodes[1][i] = arg2[i];
        nodes[2][i] = arg3[i];
        nodes[3][i] = arg4[i];
    }

    return make_polygon(nodes);
}

} // namespace geometry

} // namespace caribou

#endif //CARIBOU_GEOMETRY_QUAD_H
