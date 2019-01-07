#ifndef CARIBOU_GEOMETRY_INTERPOLATION_TRIANGLE_H
#define CARIBOU_GEOMETRY_INTERPOLATION_TRIANGLE_H

#include <Caribou/config.h>
#include <Caribou/Algebra/Vector.h>
#include <Caribou/Geometry/Interpolation/CanonicalElement.h>
#include <Caribou/Geometry/Node.h>

namespace caribou {
namespace geometry {
namespace interpolation {

/**
 * Interpolation on a quad with Lagrange polynomials of degree 1 (P1)
 *
 * v
 * ^                                                                   2
 * |                                                                   | \
 * 2                       2                    2                      9   8
 * |`\                     |`\                  | \                    |     \
 * |  `\                   |  `\                7   6                 10 (14)  7
 * |    `\                 5    `4              |     \                |         \
 * |      `\               |      `\            8  (9)  5             11 (12) (13) 6
 * |        `\             |        `\          |         \            |             \
 * 0----------1 --> u      0-----3----1         0---3---4---1          0---3---4---5---1
 *
 */
struct Triangle3 : public CanonicalElement<2, 3, Triangle3>
{
    using Index = INTEGER_TYPE;
    using Real = FLOATING_POINT_TYPE;

    static constexpr INTEGER_TYPE Dimension = 2;
    static constexpr INTEGER_TYPE NumberOfNodes = 3;

    /**
     * Compute the ith Lagrange polynomial value evaluated at local coordinates {u, v} w.r.t the triangle's interpolation node i.
     *
     * @example
     * \code{.cpp}
     * // Computes the value of node #2 Lagrange polynomial evaluated at local coordinates {-0.4, 0.2}
     * float p = Triangle3::L<2> (-0.4, 0.2);
     * \endcode
     *
     * @tparam interpolation_node_index The interpolation node id
     */
    template<INTEGER_TYPE interpolation_node_index>
    static constexpr
    Real
    L (const Real &u, const Real &v)
    {
        static_assert(interpolation_node_index >= 0 and interpolation_node_index < NumberOfNodes,
                      "The shape value can only be computed at the interpolation nodes (indices 0, 1 and 2).");

        if CONSTEXPR_IF (interpolation_node_index == 0)
            return (Real) (1 - u - v);
        else if CONSTEXPR_IF (interpolation_node_index == (Index) 1)
            return (Real) (u);
        else // interpolation_node_index == (Index) 2
            return (Real) (v);
    }

    /**
     * Compute the ith Lagrange polynomial derivatives w.r.t the local frame {dL/du, dL/dv} evaluated at local
     * coordinates {u, v} w.r.t the triangle's interpolation node i.
     *
     * @example
     * \code{.cpp}
     * // Computes the derivatives of node #2 Lagrange polynomial evaluated at local coordinates {-0.4, 0.2}
     * Vector<2, float> dp = Triangle3::dL<2> (-0.4, 0.2);
     * \endcode
     *
     * @tparam interpolation_node_index The interpolation node id
     */
    template<INTEGER_TYPE interpolation_node_index>
    static constexpr
    algebra::Vector<2, Real>
    dL (const Real & /*u*/, const Real & /*v*/)
    {
        static_assert(interpolation_node_index >= 0 and interpolation_node_index < NumberOfNodes,
                      "The shape derivatives can only be computed at the interpolation nodes (indices 0, 1 and 2).");

        if CONSTEXPR_IF (interpolation_node_index == 0)
            return {
                    -1, // dL/du
                    -1  // dL/dv
            };
        else if CONSTEXPR_IF (interpolation_node_index == (Index) 1)
            return {
                    1, // dL/du
                    0  // dL/dv
            };
        else // interpolation_node_index == (Index) 2
            return {
                    0, // dL/du
                    1  // dL/dv
            };

    }
};

} // namespace interpolation
} // namespace geometry
} // namespace caribou
#endif //CARIBOU_GEOMETRY_INTERPOLATION_TRIANGLE_H
