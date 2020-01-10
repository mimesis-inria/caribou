#ifndef CARIBOU_GEOMETRY_INTERPOLATION_QUAD_H
#define CARIBOU_GEOMETRY_INTERPOLATION_QUAD_H

#include <Caribou/config.h>
#include <Caribou/Geometry/Interpolation/CanonicalElement.h>
#include <Caribou/Geometry/Interpolation/Segment.h>
#include <Eigen/Core>

namespace caribou::geometry::interpolation {

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
struct Quad4 : public CanonicalElement<2, 4, Quad4>
{
    using Index = INTEGER_TYPE;
    using Real = FLOATING_POINT_TYPE;

    static constexpr INTEGER_TYPE Dimension = 2;
    static constexpr INTEGER_TYPE NumberOfNodes = 4;
    using LocalCoordinates = Eigen::Matrix<FLOATING_POINT_TYPE, Dimension, 1>;

    using BoundaryType = Segment2;

    static constexpr FLOATING_POINT_TYPE nodes [NumberOfNodes][Dimension] = {
    //    u,  v
        {-1, -1}, // Node 0
        {+1, -1}, // Node 1
        {+1, +1}, // Node 2
        {-1, +1}  // Node 3
    };


    static constexpr UNSIGNED_INTEGER_TYPE edges [4][2] {
        {0, 1}, // Edge 0
        {1, 2}, // Edge 1
        {2, 3}, // Edge 2
        {3, 0}  // Edge 3
    };

    static constexpr UNSIGNED_INTEGER_TYPE number_of_gauss_nodes = 4;
    static constexpr FLOATING_POINT_TYPE gauss_nodes [number_of_gauss_nodes][Dimension] {
    //          u,                w
        {-1/1.73205080757, -1/1.73205080757}, // Node 0
        {+1/1.73205080757, -1/1.73205080757}, // Node 1
        {-1/1.73205080757, +1/1.73205080757}, // Node 2
        {+1/1.73205080757, +1/1.73205080757} // Node 3
    };

    static constexpr FLOATING_POINT_TYPE gauss_weights[number_of_gauss_nodes] {1, 1, 1, 1};

    /**
     * Get the Lagrange polynomial values evaluated at local coordinates {u, v} w.r.t each quad's interpolation nodes.
     *
     * @example
     * \code{.cpp}
     * // Computes the value of node #2 Lagrange polynomial evaluated at local coordinates {-0.4, 0.2}
     * float p = Quad4::L(-0.4, 0.2)[2];
     * \endcode
     */
    static
    Eigen::Matrix<FLOATING_POINT_TYPE, NumberOfNodes, 1>
    L (const LocalCoordinates & x)
    {
        Eigen::Matrix<FLOATING_POINT_TYPE, NumberOfNodes, 1> m;
        const auto & u = x[0];
        const auto & v = x[1];

        m << (1 / 4. * (1 - u) * (1 - v)),
             (1 / 4. * (1 + u) * (1 - v)),
             (1 / 4. * (1 + u) * (1 + v)),
             (1 / 4. * (1 - u) * (1 + v));

        return m;
    }

    /**
     * Get Lagrange polynomial derivatives w.r.t the local frame {dL/du, dL/dv} evaluated at local
     * coordinates {u, v} w.r.t each quad's interpolation nodes.
     *
     * @example
     * \code{.cpp}
     * // Computes the derivatives of node #2 Lagrange polynomial evaluated at local coordinates {-0.4, 0.2}
     * Eigen::Vector2d dp = Quad4::dL (-0.4, 0.2).row(2);
     * \endcode
     */
    static
    Eigen::Matrix<FLOATING_POINT_TYPE, NumberOfNodes, Dimension, Eigen::RowMajor>
    dL (const LocalCoordinates & x)
    {
        const auto & u = x[0];
        const auto & v = x[1];

        Eigen::Matrix<FLOATING_POINT_TYPE, NumberOfNodes, Dimension, Eigen::RowMajor> m;
        //         dL/du                dL/dv
        m << -1 / 4. * (1 - v),   -1 / 4. * (1 - u),   // Node 0
             +1 / 4. * (1 - v),   -1 / 4. * (1 + u),   // Node 1
             +1 / 4. * (1 + v),   +1 / 4. * (1 + u),   // Node 2
             -1 / 4. * (1 + v),   +1 / 4. * (1 - u);   // Node 3

        return m;
    }
};

} // namespace caribou::geometry::interpolation
#endif //CARIBOU_GEOMETRY_INTERPOLATION_QUAD_H
