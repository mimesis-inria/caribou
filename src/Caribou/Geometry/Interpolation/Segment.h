#ifndef CARIBOU_GEOMETRY_INTERPOLATION_SEGMENT_H
#define CARIBOU_GEOMETRY_INTERPOLATION_SEGMENT_H

#include <Caribou/config.h>
#include <Eigen/Core>
#include <Caribou/Geometry/Interpolation/CanonicalElement.h>

namespace caribou::geometry::interpolation {

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
    static constexpr UNSIGNED_INTEGER_TYPE Dimension = 1;
    static constexpr UNSIGNED_INTEGER_TYPE NumberOfNodes = 2;
    using LocalCoordinates = Eigen::Matrix<FLOATING_POINT_TYPE, Dimension, 1>;

    template<int nRows, int nColumns, int Options=0>
    using Matrix = Eigen::Matrix<FLOATING_POINT_TYPE, nRows, nColumns, Options>;

    static constexpr FLOATING_POINT_TYPE nodes [NumberOfNodes][Dimension] = {
    //    u
        {-1}, // Node 0
        {+1}, // Node 1
    };

    static constexpr UNSIGNED_INTEGER_TYPE number_of_gauss_nodes = 1;
    static constexpr FLOATING_POINT_TYPE gauss_nodes [number_of_gauss_nodes][Dimension] {
    //   u
        {0} // Node 0
    };

    static constexpr FLOATING_POINT_TYPE gauss_weights[number_of_gauss_nodes] {2};

    /**
     * Get the Lagrange polynomial values evaluated at local coordinates {u} w.r.t each segment's interpolation nodes.
     *
     * @example
     * \code{.cpp}
     * // Computes the value of node #2 Lagrange polynomial evaluated at local coordinates {-0.4}
     * float p = Segment2::L(-0.4)[2];
     * \endcode
     */
    static
    Eigen::Matrix<FLOATING_POINT_TYPE, NumberOfNodes, 1>
    L (const LocalCoordinates & x)
    {
        Eigen::Matrix<FLOATING_POINT_TYPE, NumberOfNodes, 1> m;
        const auto  & u = x[0];

        m << 1/2. * (1 - u),
             1/2. * (1 + u);

        return m;
    }

    /**
     * Get the Lagrange polynomial derivatives w.r.t the local frame {dL/du} evaluated at local
     * coordinates {u} w.r.t each segment's interpolation nodes.
     *
     * @example
     * \code{.cpp}
     * // Computes the derivatives of node #2 Lagrange polynomial evaluated at local coordinates {-0.4}
     * float dp = Segment2::dL(-0.4)[2];
     * \endcode
     */
    static
    Matrix<NumberOfNodes, Dimension>
    dL (const LocalCoordinates & /*x*/)
    {
        Matrix<NumberOfNodes, Dimension> m;
        m << -1,
             +1;
        return m;
    }
};

struct Segment3 : public CanonicalElement<1, 3, Segment2>
{
    static constexpr UNSIGNED_INTEGER_TYPE Dimension = 1;
    static constexpr UNSIGNED_INTEGER_TYPE NumberOfNodes = 3;
    using LocalCoordinates = Eigen::Matrix<FLOATING_POINT_TYPE, Dimension, 1>;

    template<int nRows, int nColumns, int Options=0>
    using Matrix = Eigen::Matrix<FLOATING_POINT_TYPE, nRows, nColumns, Options>;

    static constexpr FLOATING_POINT_TYPE nodes [NumberOfNodes][Dimension] = {
    //    u
        {-1.0}, // Node 0
        {+1.0}, // Node 1
        {+0.5}  // Node 2
    };

    /**
     * Get the Lagrange polynomial values evaluated at local coordinates {u} w.r.t each segment's interpolation nodes.
     *
     * @example
     * \code{.cpp}
     * // Computes the value of node #2 Lagrange polynomial evaluated at local coordinates {-0.4}
     * float p = Segment3::L(-0.4)[2];
     * \endcode
     */
    static
    Eigen::Matrix<FLOATING_POINT_TYPE, NumberOfNodes, 1>
    L (const LocalCoordinates & x)
    {
        Eigen::Matrix<FLOATING_POINT_TYPE, NumberOfNodes, 1> m;
        const auto  & u = x[0];

        m << 1/2. * u * (1 - u),
             1/2. * u * (1 + u),
             1 - (u * u);

        return m;
    }

    /**
     * Get the Lagrange polynomial derivatives w.r.t the local frame {dL/du} evaluated at local
     * coordinates {u} w.r.t each segment's interpolation nodes.
     *
     * @example
     * \code{.cpp}
     * // Computes the derivatives of node #2 Lagrange polynomial evaluated at local coordinates {-0.4}
     * float dp = Segment3::dL(-0.4)[2];
     * \endcode
     */
    static
    Matrix<NumberOfNodes, Dimension>
    dL (const LocalCoordinates & x)
    {
        Matrix<NumberOfNodes, Dimension> m;
        const auto  & u = x[0];

        m << u - 1/2.,
             u + 1/2.,
             - 2 * u;
        return m;
    }
};

} // namespace caribou::geometry::interpolation
#endif //CARIBOU_GEOMETRY_INTERPOLATION_SEGMENT_H
