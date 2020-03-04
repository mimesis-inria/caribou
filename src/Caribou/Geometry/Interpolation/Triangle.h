#ifndef CARIBOU_GEOMETRY_INTERPOLATION_TRIANGLE_H
#define CARIBOU_GEOMETRY_INTERPOLATION_TRIANGLE_H

#include <Caribou/config.h>
#include <Caribou/Geometry/Interpolation/CanonicalElement.h>

namespace caribou::geometry::interpolation {

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

    static constexpr const char * name = "Triangle3";

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

    static constexpr UNSIGNED_INTEGER_TYPE number_of_gauss_nodes = 1;
    static constexpr FLOATING_POINT_TYPE gauss_nodes [number_of_gauss_nodes][Dimension] {
        //  u,    w
          {1/3., 1/3.} // Node 0
    };
    static constexpr FLOATING_POINT_TYPE gauss_weights[number_of_gauss_nodes] {1/2.};

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
    Eigen::Matrix<FLOATING_POINT_TYPE, NumberOfNodes, Dimension, Eigen::RowMajor>
    dL (const LocalCoordinates & /*x*/)
    {
        Eigen::Matrix<FLOATING_POINT_TYPE, NumberOfNodes, Dimension, Eigen::RowMajor> m;
        //    dL/du dL/dv
        m <<   -1,   -1,   // Node 0
                1,    0,   // Node 1
                0,    1;   // Node 2

        return m;
    }
};

/**
 * Interpolation on a triangle with Lagrange polynomials of degree 2 (P2)
 *
 * v
 * ^
 * |
 * 2
 * |`\
 * |  `\
 * 5    `4
 * |      `\
 * |        `\
 * 0-----3----1 --> u
 *
 */
struct Triangle6 : public CanonicalElement<2, 6, Triangle3>
{
    using Index = INTEGER_TYPE;
    using Real = FLOATING_POINT_TYPE;

    static constexpr INTEGER_TYPE Dimension = 2;
    static constexpr INTEGER_TYPE NumberOfNodes = 6;
    using LocalCoordinates = Eigen::Matrix<FLOATING_POINT_TYPE, Dimension, 1>;

    static constexpr const char * name = "Triangle6";

    static constexpr FLOATING_POINT_TYPE nodes [NumberOfNodes][Dimension] = {
    //    u,     v
        {0.0,   0.0}, // Node 0
        {1.0,   0.0}, // Node 1
        {0.0,   1.0}, // Node 2
        {0.5,   0.0}, // Node 3
        {0.5,   0.5}, // Node 4
        {0.0,   0.5}  // Node 5
    };


    static constexpr UNSIGNED_INTEGER_TYPE edges [3][3] {
        {0, 1, 3}, // Edge 0
        {1, 2, 4}, // Edge 1
        {2, 0, 5}  // Edge 2
    };

    static constexpr UNSIGNED_INTEGER_TYPE number_of_gauss_nodes = 3;
    static constexpr FLOATING_POINT_TYPE gauss_nodes [number_of_gauss_nodes][Dimension] {
    //   u,     w
        {2/3., 1/6.}, // Node 0
        {1/6., 2/3.}, // Node 1
        {1/6., 1/6.}  // Node 2
    };
    static constexpr FLOATING_POINT_TYPE gauss_weights[number_of_gauss_nodes] {
        1/6., 1/6., 1/6.
    };

    /**
     * Get the Lagrange polynomial values evaluated at local coordinates {u, v} w.r.t each triangle's interpolation nodes.
     *
     * @example
     * \code{.cpp}
     * // Computes the value of node #2 Lagrange polynomial evaluated at local coordinates {-0.4, 0.2}
     * float p = Triangle6::L({-0.4, 0.2})[2];
     * \endcode
     */
    static
    Eigen::Matrix<FLOATING_POINT_TYPE, NumberOfNodes, 1>
    L (const LocalCoordinates & x)
    {
        Eigen::Matrix<FLOATING_POINT_TYPE, NumberOfNodes, 1> m;
        const auto & u = x[0];
        const auto & v = x[1];
        const auto l = 1 - u - v;

        m << l * (2*l - 1),   // Node 0
             u * (2*u - 1),   // Node 1
             v * (2*v - 1),   // Node 2
             4 * u * l,       // Node 3
             4 * u * v,       // Node 4
             4 * v * l;       // Node 5

        return m;
    }

    /**
     * Get the Lagrange polynomial derivatives w.r.t the local frame {dL/du, dL/dv} evaluated at local
     * coordinates {u, v} w.r.t each triangle's interpolation nodes.
     *
     * @example
     * \code{.cpp}
     * // Computes the derivatives of node #2 Lagrange polynomial evaluated at local coordinates {-0.4, 0.2}
     * Eigen::Vector2d dp = Triangle6::dL (-0.4, 0.2).row(2);
     * \endcode
     */
    static
    Eigen::Matrix<FLOATING_POINT_TYPE, NumberOfNodes, Dimension, Eigen::RowMajor>
    dL (const LocalCoordinates & x)
    {
        Eigen::Matrix<FLOATING_POINT_TYPE, NumberOfNodes, Dimension, Eigen::RowMajor> m;
        const auto & u = x[0];
        const auto & v = x[1];
        const auto l = 1 - u - v;

        //       dL/du              dL/dv
        m <<   1 - 4 * l  ,       1 - 4 * l  ,   // Node 0
               4 * u - 1  ,          0       ,   // Node 1
                    0     ,       4 * v - 1  ,   // Node 2
               4 * (l - u),      -4 * u      ,   // Node 3
                  4 * v   ,       4 * u      ,   // Node 4
                - 4 * v   ,       4 * (l - v);   // Node 5

        return m;
    }
};

} // namespace caribou::geometry::interpolation
#endif //CARIBOU_GEOMETRY_INTERPOLATION_TRIANGLE_H
