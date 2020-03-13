#ifndef CARIBOU_GEOMETRY_INTERPOLATION_HEXAHEDRON_H
#define CARIBOU_GEOMETRY_INTERPOLATION_HEXAHEDRON_H

#include <Caribou/config.h>
#include <Caribou/Geometry/Interpolation/CanonicalElement.h>
#include <Caribou/Geometry/Interpolation/Quad.h>

namespace caribou::geometry::interpolation {

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

    using BoundaryType = Quad4;

    static constexpr INTEGER_TYPE Dimension = 3;
    static constexpr INTEGER_TYPE NumberOfNodes = 8;
    using LocalCoordinates = Eigen::Matrix<FLOATING_POINT_TYPE, Dimension, 1>;

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


    static constexpr UNSIGNED_INTEGER_TYPE number_of_edges = 12;
    static constexpr UNSIGNED_INTEGER_TYPE edges [number_of_edges][2] {
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

    static constexpr UNSIGNED_INTEGER_TYPE number_of_faces = 6;
    static constexpr UNSIGNED_INTEGER_TYPE faces[number_of_faces][BoundaryType::NumberOfNodes]{
            {0, 3, 2, 1}, // Face 0
            {0, 4, 7, 3}, // Face 1
            {1, 2, 6, 5}, // Face 2
            {0, 1, 5, 4}, // Face 3
            {2, 3, 7, 6}, // Face 4
            {4, 5, 6, 7}  // Face 5
    };

    static constexpr UNSIGNED_INTEGER_TYPE number_of_gauss_nodes = 8;
    static constexpr FLOATING_POINT_TYPE gauss_nodes [number_of_gauss_nodes][Dimension] {
    //                u,                 v,               w
            {-1/1.73205080757, -1/1.73205080757, -1/1.73205080757}, // Node 0
            {+1/1.73205080757, -1/1.73205080757, -1/1.73205080757}, // Node 1
            {-1/1.73205080757, +1/1.73205080757, -1/1.73205080757}, // Node 2
            {+1/1.73205080757, +1/1.73205080757, -1/1.73205080757}, // Node 3
            {-1/1.73205080757, -1/1.73205080757, +1/1.73205080757}, // Node 4
            {+1/1.73205080757, -1/1.73205080757, +1/1.73205080757}, // Node 5
            {-1/1.73205080757, +1/1.73205080757, +1/1.73205080757}, // Node 6
            {+1/1.73205080757, +1/1.73205080757, +1/1.73205080757}  // Node 7
    };

    static constexpr FLOATING_POINT_TYPE gauss_weights[number_of_gauss_nodes] {1, 1, 1, 1, 1, 1, 1, 1};

    /**
     * Get the Lagrange polynomial value evaluated at local coordinates {u, v, w} w.r.t each
     * hexa's interpolation nodes.
     *
     * @example
     * \code{.cpp}
     * // Computes the value of node #2 Lagrange polynomial evaluated at local coordinates {-0.4, 0.2, 0.1}
     * float p = Hexahedron8::L(-0.4, 0.2, 0.1)[2];
     * \endcode
     *
     * @tparam interpolation_node_index The interpolation node id
     */
    static
    Eigen::Matrix<FLOATING_POINT_TYPE, NumberOfNodes, 1>
    L (const LocalCoordinates & x)
    {
        const auto & u = x[0];
        const auto & v = x[1];
        const auto & w = x[2];

        Eigen::Matrix<FLOATING_POINT_TYPE, NumberOfNodes, 1> m;
        m << (1/8.) * (1 - u) * (1 - v) * (1 - w),
             (1/8.) * (1 + u) * (1 - v) * (1 - w),
             (1/8.) * (1 + u) * (1 + v) * (1 - w),
             (1/8.) * (1 - u) * (1 + v) * (1 - w),
             (1/8.) * (1 - u) * (1 - v) * (1 + w),
             (1/8.) * (1 + u) * (1 - v) * (1 + w),
             (1/8.) * (1 + u) * (1 + v) * (1 + w),
             (1/8.) * (1 - u) * (1 + v) * (1 + w);


        return m;
    }

    /**
     * Get the ith Lagrange polynomial derivatives w.r.t the local frame {dL/du, dL/dv, dL/dw} evaluated at local
     * coordinates {u, v, w} w.r.t each hexa's interpolation nodes.
     *
     * @example
     * \code{.cpp}
     * // Computes the derivatives of node #2 Lagrange polynomial evaluated at local coordinates {-0.4, 0.2, 0.1}
     * Eigen::Vector3d dp = Hexahedron8::dL(-0.4, 0.2, 0.1).row(2);
     * \endcode
     */
    static
    Eigen::Matrix<FLOATING_POINT_TYPE, NumberOfNodes, Dimension, Eigen::RowMajor>
    dL (const LocalCoordinates & x)
    {
        const auto & u = x[0];
        const auto & v = x[1];
        const auto & w = x[2];

        Eigen::Matrix<FLOATING_POINT_TYPE, NumberOfNodes, Dimension, Eigen::RowMajor> m;
        //            dL/du                         dL/dv                         dL/dw
        m << -1/8. * (1 - v) * (1 - w),    -1/8. * (1 - u) * (1 - w),    -1/8. * (1 - u) * (1 - v),   // Node 0
             +1/8. * (1 - v) * (1 - w),    -1/8. * (1 + u) * (1 - w),    -1/8. * (1 + u) * (1 - v),   // Node 1
             +1/8. * (1 + v) * (1 - w),    +1/8. * (1 + u) * (1 - w),    -1/8. * (1 + u) * (1 + v),   // Node 2
             -1/8. * (1 + v) * (1 - w),    +1/8. * (1 - u) * (1 - w),    -1/8. * (1 - u) * (1 + v),   // Node 3
             -1/8. * (1 - v) * (1 + w),    -1/8. * (1 - u) * (1 + w),    +1/8. * (1 - u) * (1 - v),   // Node 4
             +1/8. * (1 - v) * (1 + w),    -1/8. * (1 + u) * (1 + w),    +1/8. * (1 + u) * (1 - v),   // Node 5
             +1/8. * (1 + v) * (1 + w),    +1/8. * (1 + u) * (1 + w),    +1/8. * (1 + u) * (1 + v),   // Node 6
             -1/8. * (1 + v) * (1 + w),    +1/8. * (1 - u) * (1 + w),    +1/8. * (1 - u) * (1 + v);   // Node 7
        return m;
    }
};

} // namespace caribou::geometry::interpolation
#endif //CARIBOU_GEOMETRY_INTERPOLATION_HEXAHEDRON_H
