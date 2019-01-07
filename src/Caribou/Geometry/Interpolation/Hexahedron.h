#ifndef CARIBOU_GEOMETRY_INTERPOLATION_HEXAHEDRON_H
#define CARIBOU_GEOMETRY_INTERPOLATION_HEXAHEDRON_H

#include <Caribou/config.h>
#include <Caribou/Algebra/Vector.h>
#include <Caribou/Geometry/Interpolation/CanonicalElement.h>
#include <Caribou/Geometry/Node.h>

namespace caribou {
namespace geometry {
namespace interpolation {

/**
 * Interpolation on an hexahedron with Lagrange polynomials of degree 1 (P1)
 *
 * Hexahedron:             Hexahedron20:          Hexahedron27:
 *
 *        v
 * 3----------2            3----13----2           3----13----2
 * |\     ^   |\           |\         |\          |\         |\
 * | \    |   | \          | 15       | 14        |15    24  | 14
 * |  \   |   |  \         9  \       11 \        9  \ 20    11 \
 * |   7------+---6        |   7----19+---6       |   7----19+---6
 * |   |  +-- |-- | -> u   |   |      |   |       |22 |  26  | 23|
 * 0---+---\--1   |        0---+-8----1   |       0---+-8----1   |
 *  \  |    \  \  |         \  17      \  18       \ 17    25 \  18
 *   \ |     \  \ |         10 |        12|        10 |  21    12|
 *    \|      w  \|           \|         \|          \|         \|
 *     4----------5            4----16----5           4----16----5
 *
 */
struct Hexahedron8 : public CanonicalElement<3, 8, Quad4>
{
    using Index = INTEGER_TYPE;
    using Real = FLOATING_POINT_TYPE;

    static constexpr INTEGER_TYPE Dimension = 3;
    static constexpr INTEGER_TYPE NumberOfNodes = 8;

    static constexpr std::array<algebra::Vector<3, Real>, 8> nodes {{
    //    u,  v,  w
        {-1, -1, -1}, // Node 0
        {+1, -1, -1}, // Node 1
        {+1, +1, -1}, // Node 2
        {-1, +1, -1}, // Node 3
        {-1, -1, +1}, // Node 4
        {+1, -1, +1}, // Node 5
        {+1, +1, +1}, // Node 6
        {-1, +1, +1}  // Node 7
    }};

    /**
     * Compute the ith Lagrange polynomial value evaluated at local coordinates {u, v, w} w.r.t the
     * hexa's interpolation node i.
     *
     * @example
     * \code{.cpp}
     * // Computes the value of node #2 Lagrange polynomial evaluated at local coordinates {-0.4, 0.2, 0.1}
     * float p = Hexahedron8::L<2> (-0.4, 0.2, 0.1);
     * \endcode
     *
     * @tparam interpolation_node_index The interpolation node id
     */
    template<INTEGER_TYPE interpolation_node_index>
    static constexpr
    Real
    L (const Real &u, const Real &v, const Real & w)
    {
        static_assert(interpolation_node_index >= 0 and interpolation_node_index < NumberOfNodes,
                      "The shape value can only be computed at the interpolation nodes.");

        constexpr Real ui = nodes[interpolation_node_index][0];
        constexpr Real vi = nodes[interpolation_node_index][1];
        constexpr Real wi = nodes[interpolation_node_index][2];

        return (Real) (1/8.) * (1 + ui*u) * (1 + vi*v) * (1 + wi*w);
    }

    /**
     * Compute the ith Lagrange polynomial derivatives w.r.t the local frame {dL/du, dL/dv, dL/dw} evaluated at local
     * coordinates {u, v, w} w.r.t the hexa's interpolation node i.
     *
     * @example
     * \code{.cpp}
     * // Computes the derivatives of node #2 Lagrange polynomial evaluated at local coordinates {-0.4, 0.2, 0.1}
     * Vector<3, float> dp = Hexahedron8::dL<2> (-0.4, 0.2, 0.1);
     * \endcode
     *
     * @tparam interpolation_node_index The interpolation node id
     */
    template<INTEGER_TYPE interpolation_node_index>
    static constexpr
    algebra::Vector<3, Real>
    dL (const Real &u, const Real &v, const Real & w)
    {
        static_assert(interpolation_node_index >= 0 and interpolation_node_index < NumberOfNodes,
                      "The shape derivatives can only be computed at the interpolation nodes.");

        constexpr Real ui = nodes[interpolation_node_index][0];
        constexpr Real vi = nodes[interpolation_node_index][1];
        constexpr Real wi = nodes[interpolation_node_index][2];

        return {
                (Real)1 / 8 *     1      * (1 + vi*v) * (1 + wi*w),
                (Real)1 / 8 * (1 + ui*u) *     1      * (1 + wi*w),
                (Real)1 / 8 * (1 + ui*u) * (1 + vi*v) *     1
        };
    }
};

} // namespace interpolation
} // namespace geometry
} // namespace caribou
#endif //CARIBOU_GEOMETRY_INTERPOLATION_HEXAHEDRON_H
