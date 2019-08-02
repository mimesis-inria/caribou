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
 * Interpolation on a triangle with Lagrange polynomials of degree 1 (P1)
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
    using LocalCoordinates = Eigen::Matrix<FLOATING_POINT_TYPE, Dimension, 1>;

    static constexpr FLOATING_POINT_TYPE nodes [NumberOfNodes][Dimension] = {
    //   u, v
        {0, 0}, // Node 0
        {1, 0}, // Node 1
        {0, 1}  // Node 2
    };


    static constexpr UNSIGNED_INTEGER_TYPE edges [3][2] {
        {0, 1}, // Edge 0
        {1, 2}, // Edge 1
        {2, 0}  // Edge 2
    };

    /**
     * Get the Lagrange polynomial values evaluated at local coordinates {u, v} w.r.t each triangle's interpolation nodes.
     *
     * @example
     * \code{.cpp}
     * // Computes the value of node #2 Lagrange polynomial evaluated at local coordinates {-0.4, 0.2}
     * float p = Triangle3::L({-0.4, 0.2})[2];
     * \endcode
     */
    static
    Eigen::Matrix<FLOATING_POINT_TYPE, NumberOfNodes, 1>
    L (const LocalCoordinates & x)
    {
        Eigen::Matrix<FLOATING_POINT_TYPE, NumberOfNodes, 1> m;
        const auto & u = x[0];
        const auto & v = x[1];

        m << 1 - u - v,
             u,
             v;

        return m;
    }

    /**
     * Get the Lagrange polynomial derivatives w.r.t the local frame {dL/du, dL/dv} evaluated at local
     * coordinates {u, v} w.r.t each triangle's interpolation nodes.
     *
     * @example
     * \code{.cpp}
     * // Computes the derivatives of node #2 Lagrange polynomial evaluated at local coordinates {-0.4, 0.2}
     * Eigen::Vector2d dp = Triangle3::dL (-0.4, 0.2).row(2);
     * \endcode
     */
    static
    Eigen::Matrix<FLOATING_POINT_TYPE, NumberOfNodes, Dimension>
    dL (const LocalCoordinates & /*x*/)
    {
        Eigen::Matrix<FLOATING_POINT_TYPE, NumberOfNodes, Dimension> m;
        //    dL/du dL/dv
        m <<   -1,   -1,   // Node 0
                1,    0,   // Node 1
                0,    1;   // Node 2

        return m;
    }
};

} // namespace interpolation
} // namespace geometry
} // namespace caribou
#endif //CARIBOU_GEOMETRY_INTERPOLATION_TRIANGLE_H
