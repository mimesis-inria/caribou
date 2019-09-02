#ifndef CARIBOU_GEOMETRY_INTERPOLATION_TETRAHEDRON_H
#define CARIBOU_GEOMETRY_INTERPOLATION_TETRAHEDRON_H

#include <Caribou/config.h>
#include <Caribou/Geometry/Interpolation/CanonicalElement.h>
#include <Caribou/Geometry/Interpolation/Triangle.h>

namespace caribou::geometry::interpolation {

/**
 * Interpolation on a tetrahedron with Lagrange polynomials of degree 1 (P1)
 *
 *
 *                    v
 *                  .
 *                ,/
 *               /
 *            2
 *          ,/|`\
 *        ,/  |  `\
 *      ,/    '.   `\
 *    ,/       |     `\
 *  ,/         |       `\
 * 0-----------'.--------1 --> u
 *  `\.         |      ,/
 *     `\.      |    ,/
 *        `\.   '. ,/
 *           `\. |/
 *              `3
 *                 `\.
 *                    ` w
 *
 */
struct Tetrahedron4 : public CanonicalElement<3, 4, Tetrahedron4>
{
    using Index = INTEGER_TYPE;
    using Real = FLOATING_POINT_TYPE;

    using FaceType = Triangle3;

    static constexpr INTEGER_TYPE Dimension = 3;
    static constexpr INTEGER_TYPE NumberOfNodes = 4;
    using LocalCoordinates = Eigen::Matrix<FLOATING_POINT_TYPE, Dimension, 1>;

    static constexpr const char * name = "Tetrahedron4";

    static constexpr FLOATING_POINT_TYPE nodes [NumberOfNodes][Dimension] = {
        //   u, v, w
            {0, 0, 0}, // Node 0
            {1, 0, 0}, // Node 1
            {0, 1, 0}, // Node 2
            {0, 0, 1}  // Node 3
    };


    static constexpr UNSIGNED_INTEGER_TYPE edges [6][2] {
        {0, 1}, // Edge 0
        {1, 2}, // Edge 1
        {2, 0}, // Edge 2
        {3, 0}, // Edge 3
        {3, 2}, // Edge 4
        {3, 1}  // Edge 5
    };

    static constexpr UNSIGNED_INTEGER_TYPE faces[4][3]{
        {0, 2, 1}, // Face 0
        {0, 1, 3}, // Face 1
        {0, 3, 2}, // Face 2
        {3, 1, 2}  // Face 3
    };

    static constexpr UNSIGNED_INTEGER_TYPE number_of_gauss_nodes = 1;
    static constexpr FLOATING_POINT_TYPE gauss_nodes [number_of_gauss_nodes][Dimension] {
        //    u,    v,    w
            {1/4., 1/4., 1/4.}, // Node 0
    };

    static constexpr FLOATING_POINT_TYPE gauss_weights[number_of_gauss_nodes] {1/6.};

    /**
     * Get the Lagrange polynomial value evaluated at local coordinates {u, v, w} w.r.t each
     * tetra's interpolation nodes.
     *
     * @example
     * \code{.cpp}
     * // Computes the value of node #2 Lagrange polynomial evaluated at local coordinates {-0.4, 0.2, 0.1}
     * float p = Tetrahedron4::L(-0.4, 0.2, 0.1)[2];
     * \endcode
     */
    static
    Eigen::Matrix<FLOATING_POINT_TYPE, NumberOfNodes, 1>
    L (const LocalCoordinates & x)
    {
        const auto & u = x[0];
        const auto & v = x[1];
        const auto & w = x[2];

        Eigen::Matrix<FLOATING_POINT_TYPE, NumberOfNodes, 1> m;
        m << 1 - u -v - w,
            u,
            v,
            w;

        return m;
    }

    /**
     * Get the ith Lagrange polynomial derivatives w.r.t the local frame {dL/du, dL/dv, dL/dw} evaluated at local
     * coordinates {u, v, w} w.r.t each tetra's interpolation nodes.
     *
     * @example
     * \code{.cpp}
     * // Computes the derivatives of node #2 Lagrange polynomial evaluated at local coordinates {-0.4, 0.2, 0.1}
     * Eigen::Vector3d dp = Tetrahedron4::dL(-0.4, 0.2, 0.1).row(2);
     * \endcode
     */
    static
    Eigen::Matrix<FLOATING_POINT_TYPE, NumberOfNodes, Dimension, Eigen::RowMajor>
    dL (const LocalCoordinates & /* x */)
    {
        Eigen::Matrix<FLOATING_POINT_TYPE, NumberOfNodes, Dimension, Eigen::RowMajor> m;
        //  dL/du  dL/dv  dL/dw
        m << -1,    -1,    -1,   // Node 0
              1,     0,     0,   // Node 1
              0,     1,     0,   // Node 1
              0,     0,     1;   // Node 3
        return m;
    }
};

/**
 * Interpolation on a tetrahedron with Lagrange polynomials of degree 2 (P2)
 *
 *
 *                    v
 *                  .
 *                ,/
 *               /
 *            2
 *          ,/|`\
 *        ,/  |  `\
 *      ,6    '.   `5
 *    ,/       8     `\
 *  ,/         |       `\
 * 0--------4--'.--------1 --> u
 *  `\.         |      ,/
 *     `\.      |    ,9
 *        `7.   '. ,/
 *           `\. |/
 *              `3
 *                 `\.
 *                    ` w
 *
 */
struct Tetrahedron10 : public CanonicalElement<3, 10, Tetrahedron4>
{
    using Index = INTEGER_TYPE;
    using Real = FLOATING_POINT_TYPE;

    using FaceType = Triangle3;

    static constexpr INTEGER_TYPE Dimension = 3;
    static constexpr INTEGER_TYPE NumberOfNodes = 10;
    using LocalCoordinates = Eigen::Matrix<FLOATING_POINT_TYPE, Dimension, 1>;

    static constexpr const char * name = "Tetrahedron10";

    static constexpr FLOATING_POINT_TYPE nodes [NumberOfNodes][Dimension] = {
        //   u, v, w
        {0.0, 0.0, 0.0}, // Node 0
        {1.0, 0.0, 0.0}, // Node 1
        {0.0, 1.0, 0.0}, // Node 2
        {0.0, 0.0, 1.0}, // Node 3
        {0.5, 0.0, 0.0}, // Node 4
        {0.5, 0.5, 0.0}, // Node 5
        {0.0, 0.5, 0.0}, // Node 6
        {0.0, 0.0, 0.5}, // Node 7
        {0.0, 0.5, 0.5}, // Node 8
        {0.5, 0.0, 0.5}  // Node 9
    };


    static constexpr UNSIGNED_INTEGER_TYPE edges [6][3] {
        {0, 1, 4}, // Edge 0
        {1, 2, 5}, // Edge 1
        {2, 0, 6}, // Edge 2
        {3, 0, 7}, // Edge 3
        {3, 2, 8}, // Edge 4
        {3, 1, 9}  // Edge 5
    };

    static constexpr UNSIGNED_INTEGER_TYPE faces[4][6]{
        {0, 2, 1, 6, 5, 4}, // Face 0
        {0, 1, 3, 4, 9, 7}, // Face 1
        {0, 3, 2, 7, 8, 6}, // Face 2
        {3, 1, 2, 9, 5, 8}  // Face 3
    };

    static constexpr UNSIGNED_INTEGER_TYPE number_of_gauss_nodes = 4;
    static constexpr FLOATING_POINT_TYPE gauss_nodes [number_of_gauss_nodes][Dimension] {
    //           u,                  v,                   w
        {0.1381966011250105, 0.1381966011250105, 0.1381966011250105}, // Node 0
        {0.1381966011250105, 0.1381966011250105, 0.5854101966249685}, // Node 1
        {0.1381966011250105, 0.5854101966249685, 0.1381966011250105}, // Node 2
        {0.5854101966249685, 0.1381966011250105, 0.1381966011250105}, // Node 3
    };

    static constexpr FLOATING_POINT_TYPE gauss_weights[number_of_gauss_nodes] {
        0.25 / 6,
        0.25 / 6,
        0.25 / 6,
        0.25 / 6
    };

    /**
     * Get the Lagrange polynomial value evaluated at local coordinates {u, v, w} w.r.t each
     * tetra's interpolation nodes.
     *
     * @example
     * \code{.cpp}
     * // Computes the value of node #2 Lagrange polynomial evaluated at local coordinates {-0.4, 0.2, 0.1}
     * float p = Tetrahedron10::L(-0.4, 0.2, 0.1)[2];
     * \endcode
     */
    static
    Eigen::Matrix<FLOATING_POINT_TYPE, NumberOfNodes, 1>
    L (const LocalCoordinates & x)
    {
        const auto & u = x[0];
        const auto & v = x[1];
        const auto & w = x[2];
        const auto   l = 1 - u - v - w;

        Eigen::Matrix<FLOATING_POINT_TYPE, NumberOfNodes, 1> m;
        m << l * (2*l - 1),   // Node 0
             u * (2*u - 1),   // Node 1
             v * (2*v - 1),   // Node 2
             w * (2*w - 1),   // Node 3
             4 * l * u,       // Node 4
             4 * u * v,       // Node 5
             4 * v * l,       // Node 6
             4 * l * w,       // Node 7
             4 * u * w,       // Node 8
             4 * v * w;       // Node 9

        return m;
    }

    /**
     * Get the ith Lagrange polynomial derivatives w.r.t the local frame {dL/du, dL/dv, dL/dw} evaluated at local
     * coordinates {u, v, w} w.r.t each tetra's interpolation nodes.
     *
     * @example
     * \code{.cpp}
     * // Computes the derivatives of node #2 Lagrange polynomial evaluated at local coordinates {-0.4, 0.2, 0.1}
     * Eigen::Vector3d dp = Tetrahedron10::dL(-0.4, 0.2, 0.1).row(2);
     * \endcode
     */
    static
    Eigen::Matrix<FLOATING_POINT_TYPE, NumberOfNodes, Dimension, Eigen::RowMajor>
    dL (const LocalCoordinates &  x)
    {
        const auto & u = x[0];
        const auto & v = x[1];
        const auto & w = x[2];
        const auto   l = 1 - u - v - w;

        Eigen::Matrix<FLOATING_POINT_TYPE, NumberOfNodes, Dimension, Eigen::RowMajor> m;
        //       dL/du              dL/dv              dL/dw
        m <<   1 - (4 * l),       1 - (4 * l),       1 - (4 * l),   // Node 0
              (4 * u) - 1 ,          0       ,          0       ,   // Node 1
                    0     ,      (4 * v) - 1 ,          0       ,   // Node 2
                    0     ,          0       ,      (4 * w) - 1 ,   // Node 3
               4 * (l - u),      -4 * u      ,      -4 * u      ,   // Node 4
                4 * v     ,       4 * u      ,          0       ,   // Node 5
               -4 * v     ,       4 * (l - v),      -4 * v      ,   // Node 6
               -4 * w     ,      -4 * w      ,       4 * (l - w),   // Node 7
                4 * w     ,          0       ,       4 * u      ,   // Node 8
                    0     ,       4 * w      ,       4 * v      ;   // Node 9
        return m;
    }
};

} // namespace caribou::geometry::interpolation
#endif //CARIBOU_GEOMETRY_INTERPOLATION_TETRAHEDRON_H
