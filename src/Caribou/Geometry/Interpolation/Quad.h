#ifndef CARIBOU_GEOMETRY_INTERPOLATION_QUAD_H
#define CARIBOU_GEOMETRY_INTERPOLATION_QUAD_H

#include <Caribou/config.h>
#include <Caribou/Algebra/Vector.h>
#include <Caribou/Geometry/Interpolation/InterpolationElement.h>
#include <Caribou/Geometry/Node.h>

namespace caribou {
namespace geometry {
namespace interpolation {

/**
 * Interpolation on a quad with Lagrange polynomials of degree 1 (P1)
 *
 *       v
 *       ^
 *       |
 * 3-----------2          3-----6-----2           3-----6-----2
 * |     |     |          |           |           |           |
 * |     |     |          |           |           |           |
 * |     +---- | --> u    7           5           7     8     5
 * |           |          |           |           |           |
 * |           |          |           |           |           |
 * 0-----------1          0-----4-----1           0-----4-----1
 */
struct Quad4
{
    using Index = INTEGER_TYPE;
    using Real = FLOATING_POINT_TYPE;

    static constexpr INTEGER_TYPE Dimension = 2;
    static constexpr INTEGER_TYPE NumberOfNodes = 4;

    /**
     * Compute the ith Lagrange polynomial value evaluated at local coordinates {u, v} w.r.t the quad's interpolation node i.
     *
     * @example
     * \code{.cpp}
     * // Computes the value of node #2 Lagrange polynomial evaluated at local coordinates {-0.4, 0.2}
     * float p = Quad4::L<2> (-0.4, 0.2);
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
                      "The shape value can only be computed at the interpolation nodes (indices 0, 1, 2 and 3).");

        if CONSTEXPR_IF (interpolation_node_index == 0)
            return (Real) (1 / 4. * (1 - u) * (1 - v));
        else if CONSTEXPR_IF (interpolation_node_index == (Index) 1)
            return (Real) (1 / 4. * (1 + u) * (1 - v));
        else if CONSTEXPR_IF (interpolation_node_index == (Index) 2)
            return (Real) (1 / 4. * (1 + u) * (1 + v));
        else // interpolation_node_index == (Index) 3
            return (Real) (1 / 4. * (1 - u) * (1 + v));
    }

    /**
     * Compute the ith Lagrange polynomial derivatives w.r.t the local frame {dL/du, dL/dv} evaluated at local
     * coordinates {u, v} w.r.t the quad's interpolation node i.
     *
     * @example
     * \code{.cpp}
     * // Computes the derivatives of node #2 Lagrange polynomial evaluated at local coordinates {-0.4, 0.2}
     * Vector<2, float> dp = Quad4::dL<2> (-0.4, 0.2);
     * \endcode
     *
     * @tparam interpolation_node_index The interpolation node id
     */
    template<INTEGER_TYPE interpolation_node_index>
    static constexpr
    algebra::Vector<2, Real>
    dL (const Real &u, const Real &v)
    {
        static_assert(interpolation_node_index >= 0 and interpolation_node_index < NumberOfNodes,
                      "The shape derivatives can only be computed at the interpolation nodes (indices 0, 1, 2 and 3).");

        if CONSTEXPR_IF (interpolation_node_index == 0)
            return {
                    -1 / 4. * (1 - v),
                    -1 / 4. * (1 - u)
            };
        else if CONSTEXPR_IF (interpolation_node_index == (Index) 1)
            return {
                    +1 / 4. * (1 - v),
                    -1 / 4. * (1 + u)
            };
        else if CONSTEXPR_IF (interpolation_node_index == (Index) 2)
            return {
                    +1 / 4. * (1 + v),
                    +1 / 4. * (1 + u)
            };
        else // interpolation_node_index == (Index) 3
            return {
                    -1 / 4. * (1 + v),
                    +1 / 4. * (1 - u)
            };
    }
};

} // namespace interpolation
} // namespace geometry
} // namespace caribou
#endif //CARIBOU_GEOMETRY_INTERPOLATION_QUAD_H
