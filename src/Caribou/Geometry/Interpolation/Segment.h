#ifndef CARIBOU_GEOMETRY_INTERPOLATION_SEGMENT_H
#define CARIBOU_GEOMETRY_INTERPOLATION_SEGMENT_H

#include <Caribou/config.h>
#include <Caribou/Algebra/Vector.h>
#include <Caribou/Geometry/Interpolation/CanonicalElement.h>
#include <Caribou/Geometry/Node.h>

namespace caribou {
namespace geometry {
namespace interpolation {

/**
 * Interpolation on a segment with Lagrange polynomials of degree 1 (P1)
 *
 *
 * P1 : 0-----+-----1 --> u
 * P2 : 0-----2-----1 --> u
 * P3 : 0--2--+--3--1--> u
 *
 *
 */
struct Segment2 : public CanonicalElement<1, 2, Segment2>
{
    using Index = INTEGER_TYPE;
    using Real = FLOATING_POINT_TYPE;

    static constexpr INTEGER_TYPE Dimension = 1;
    static constexpr INTEGER_TYPE NumberOfNodes = 2;

    /**
     * Compute the ith Lagrange polynomial value evaluated at local coordinates {u} w.r.t the segment's interpolation node i.
     *
     * @example
     * \code{.cpp}
     * // Computes the value of node #2 Lagrange polynomial evaluated at local coordinates {-0.4}
     * float p = Segment2::L<2> (-0.4);
     * \endcode
     *
     * @tparam interpolation_node_index The interpolation node id
     */
    template<INTEGER_TYPE interpolation_node_index>
    static constexpr
    Real
    L (const Real &u)
    {
        static_assert(interpolation_node_index >= 0 and interpolation_node_index < NumberOfNodes,
                      "The shape value can only be computed at the interpolation nodes (indices 0 and 1).");

        if CONSTEXPR_IF (interpolation_node_index == 0)
            return (Real) 1/2. * (1 - u);
        else // interpolation_node_index == (Index) 1
            return (Real) 1/2. * (1 + u);
    }

    /**
     * Compute the ith Lagrange polynomial derivatives w.r.t the local frame {dL/du} evaluated at local
     * coordinates {u} w.r.t the segment's interpolation node i.
     *
     * @example
     * \code{.cpp}
     * // Computes the derivatives of node #2 Lagrange polynomial evaluated at local coordinates {-0.4}
     * float dp = Segment2::dL<2> (-0.4);
     * \endcode
     *
     * @tparam interpolation_node_index The interpolation node id
     */
    template<INTEGER_TYPE interpolation_node_index>
    static constexpr
    algebra::Vector<1, Real>
    dL (const Real & /*u*/)
    {
        static_assert(interpolation_node_index >= 0 and interpolation_node_index < NumberOfNodes,
                      "The shape derivatives can only be computed at the interpolation nodes (indices 0 and 1).");

        if CONSTEXPR_IF (interpolation_node_index == 0)
            return {
                    -1 // dL/du
            };
        else // interpolation_node_index == (Index) 1
            return {
                    1 // dL/du
            };

    }
};

} // namespace interpolation
} // namespace geometry
} // namespace caribou
#endif //CARIBOU_GEOMETRY_INTERPOLATION_SEGMENT_H
