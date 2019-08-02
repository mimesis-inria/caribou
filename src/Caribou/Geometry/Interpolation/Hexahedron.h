#ifndef CARIBOU_GEOMETRY_INTERPOLATION_HEXAHEDRON_H
#define CARIBOU_GEOMETRY_INTERPOLATION_HEXAHEDRON_H

#include <Caribou/config.h>
#include <Caribou/Algebra/Vector.h>
#include <Caribou/Geometry/Interpolation/CanonicalElement.h>
#include <Caribou/Geometry/Interpolation/Quad.h>
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
struct Hexahedron8 : public CanonicalElement<3, 8, Hexahedron8>
{
    using Index = INTEGER_TYPE;
    using Real = FLOATING_POINT_TYPE;

    using QuadType = Quad4;

    static constexpr INTEGER_TYPE Dimension = 3;
    static constexpr INTEGER_TYPE NumberOfNodes = 8;

    static constexpr FLOATING_POINT_TYPE nodes [NumberOfNodes][Dimension] = {
        //    u,  v,  w
        {-1, -1, -1}, // Node 0
        {+1, -1, -1}, // Node 1
        {+1, +1, -1}, // Node 2
        {-1, +1, -1}, // Node 3
        {-1, -1, +1}, // Node 4
        {+1, -1, +1}, // Node 5
        {+1, +1, +1}, // Node 6
        {-1, +1, +1}  // Node 7
    };


    static constexpr UNSIGNED_INTEGER_TYPE edges [12][2] {
        {0, 1}, // Edge 0
        {1, 2}, // Edge 1
        {2, 3}, // Edge 2
        {3, 0}, // Edge 3
        {0, 4}, // Edge 4
        {3, 7}, // Edge 5
        {2, 6}, // Edge 6
        {1, 5}, // Edge 7
        {4, 5}, // Edge 8
        {5, 6}, // Edge 9
        {6, 7}, // Edge 10
        {7, 4}  // Edge 11
    };

    static constexpr UNSIGNED_INTEGER_TYPE faces [6][4] {
        {0, 3, 2, 1}, // Face 0
        {0, 4, 7, 3}, // Face 1
        {1, 2, 6, 5}, // Face 2
        {0, 1, 5, 4}, // Face 3
        {2, 3, 7, 6}, // Face 4
        {4, 5, 6, 7}  // Face 5
    };

    static constexpr std::array<caribou::geometry::Node<3>, 8> gauss_nodes {{
    //    u,  v,  w
        {-1/1.73205080757, -1/1.73205080757, -1/1.73205080757}, // Node 0
        {+1/1.73205080757, -1/1.73205080757, -1/1.73205080757}, // Node 1
        {+1/1.73205080757, +1/1.73205080757, -1/1.73205080757}, // Node 2
        {-1/1.73205080757, +1/1.73205080757, -1/1.73205080757}, // Node 3
        {-1/1.73205080757, -1/1.73205080757, +1/1.73205080757}, // Node 4
        {+1/1.73205080757, -1/1.73205080757, +1/1.73205080757}, // Node 5
        {+1/1.73205080757, +1/1.73205080757, +1/1.73205080757}, // Node 6
        {-1/1.73205080757, +1/1.73205080757, +1/1.73205080757}  // Node 7
    }};

    static constexpr std::array<Real, 8> gauss_weights {{1, 1, 1, 1, 1, 1, 1, 1}};

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
                (Real)1 / 8 *     ui     * (1 + vi*v) * (1 + wi*w),
                (Real)1 / 8 * (1 + ui*u) *     vi     * (1 + wi*w),
                (Real)1 / 8 * (1 + ui*u) * (1 + vi*v) *     wi
        };
    }
};

} // namespace interpolation
} // namespace geometry
} // namespace caribou
#endif //CARIBOU_GEOMETRY_INTERPOLATION_HEXAHEDRON_H
