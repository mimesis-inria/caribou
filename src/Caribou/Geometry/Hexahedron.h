#pragma once

#include <Caribou/config.h>
#include <Caribou/Geometry/BaseHexahedron.h>
#include <Caribou/Geometry/Quad.h>
#include <Eigen/Core>

namespace caribou::geometry {


template<>
struct traits<Hexahedron <Linear>> {
    static constexpr UNSIGNED_INTEGER_TYPE CanonicalDimension = 3;
    static constexpr UNSIGNED_INTEGER_TYPE Dimension = 3;
    static constexpr UNSIGNED_INTEGER_TYPE NumberOfNodesAtCompileTime = 8;
    static constexpr UNSIGNED_INTEGER_TYPE NumberOfGaussNodesAtCompileTime = 8;

    using BoundaryElementType = Quad<3, Linear>;
    static constexpr UNSIGNED_INTEGER_TYPE NumberOfBoundaryElementsAtCompileTime = 6;
};

/**
 * Linear Hexahedron
 *
 *        v
 * 3----------2
 * |\     ^   |\
 * | \    |   | \
 * |  \   |   |  \
 * |   7------+---6
 * |   |  +-- |-- | -> u
 * 0---+---\--1   |
 *  \  |    \  \  |
 *   \ |     \  \ |
 *    \|      w  \|
 *     4----------5
 *
 */
template<>
struct Hexahedron <Linear> : public BaseHexahedron<Hexahedron <Linear>> {
    // Types
    using Base = BaseHexahedron<Hexahedron <Linear>>;
    using LocalCoordinates = typename Base::LocalCoordinates;
    using WorldCoordinates = typename Base::WorldCoordinates;

    using GaussNode = typename Base::GaussNode;

    template <UNSIGNED_INTEGER_TYPE Dim>
    using Vector = typename Base::template Vector<Dim>;

    template <UNSIGNED_INTEGER_TYPE Rows, UNSIGNED_INTEGER_TYPE Cols>
    using Matrix = typename Base::template Matrix<Rows, Cols>;

    // Constants
    static constexpr auto CanonicalDimension = Base::CanonicalDimension;
    static constexpr auto Dimension = Base::Dimension;
    static constexpr auto NumberOfNodesAtCompileTime = Base::NumberOfNodesAtCompileTime;
    static constexpr auto NumberOfGaussNodesAtCompileTime = Base::NumberOfGaussNodesAtCompileTime;

    // Constructors
    using Base::Base;
    Hexahedron() : Base() {
        this->p_nodes.row(0) = WorldCoordinates(-1, -1, -1);
        this->p_nodes.row(1) = WorldCoordinates(+1, -1, -1);
        this->p_nodes.row(2) = WorldCoordinates(+1, +1, -1);
        this->p_nodes.row(3) = WorldCoordinates(-1, +1, -1);
        this->p_nodes.row(4) = WorldCoordinates(-1, -1, +1);
        this->p_nodes.row(5) = WorldCoordinates(+1, -1, +1);
        this->p_nodes.row(6) = WorldCoordinates(+1, +1, +1);
        this->p_nodes.row(7) = WorldCoordinates(-1, +1, +1);
    };


private:
    // Implementations
    friend struct Element<Hexahedron <Linear>>;
    friend struct BaseHexahedron<Hexahedron <Linear>>;
    inline auto get_L(const LocalCoordinates & xi) const {
        const auto  & u = xi[0];
        const auto  & v = xi[1];
        const auto  & w = xi[2];

        Vector<NumberOfNodesAtCompileTime> m;
        m << (1/8.) * (1 - u) * (1 - v) * (1 - w),
             (1/8.) * (1 + u) * (1 - v) * (1 - w),
             (1/8.) * (1 + u) * (1 + v) * (1 - w),
             (1/8.) * (1 - u) * (1 + v) * (1 - w),
             (1/8.) * (1 - u) * (1 - v) * (1 + w),
             (1/8.) * (1 + u) * (1 - v) * (1 + w),
             (1/8.) * (1 + u) * (1 + v) * (1 + w),
             (1/8.) * (1 - u) * (1 + v) * (1 + w);
        return m;
    };

    inline auto get_dL(const LocalCoordinates & xi) const {
        const auto  & u = xi[0];
        const auto  & v = xi[1];
        const auto  & w = xi[2];

        Matrix<NumberOfNodesAtCompileTime, CanonicalDimension> m;
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
    };

    inline auto get_gauss_nodes() const -> const auto & {
        using Weight = FLOATING_POINT_TYPE;
        static const std::vector<GaussNode> gauss_nodes {
            GaussNode {LocalCoordinates(-1/sqrt(3), -1/sqrt(3), -1/sqrt(3)), Weight(1)}, // Node 0
            GaussNode {LocalCoordinates(+1/sqrt(3), -1/sqrt(3), -1/sqrt(3)), Weight(1)}, // Node 1
            GaussNode {LocalCoordinates(-1/sqrt(3), +1/sqrt(3), -1/sqrt(3)), Weight(1)}, // Node 2
            GaussNode {LocalCoordinates(+1/sqrt(3), +1/sqrt(3), -1/sqrt(3)), Weight(1)}, // Node 3
            GaussNode {LocalCoordinates(-1/sqrt(3), -1/sqrt(3), +1/sqrt(3)), Weight(1)}, // Node 4
            GaussNode {LocalCoordinates(+1/sqrt(3), -1/sqrt(3), +1/sqrt(3)), Weight(1)}, // Node 5
            GaussNode {LocalCoordinates(-1/sqrt(3), +1/sqrt(3), +1/sqrt(3)), Weight(1)}, // Node 6
            GaussNode {LocalCoordinates(+1/sqrt(3), +1/sqrt(3), +1/sqrt(3)), Weight(1)}  // Node 7
        };
        return gauss_nodes;
    }

    inline auto get_boundary_elements_nodes() const -> const auto & {
        static const std::array<std::array<UNSIGNED_INTEGER_TYPE, 4>, 6> indices = {{
            {0, 3, 2, 1}, // Face 0
            {0, 4, 7, 3}, // Face 1
            {1, 2, 6, 5}, // Face 2
            {0, 1, 5, 4}, // Face 3
            {2, 3, 7, 6}, // Face 4
            {4, 5, 6, 7}  // Face 5
        }};
        return indices;
    }
};


template<>
struct traits<Hexahedron <Quadratic>> {
    static constexpr UNSIGNED_INTEGER_TYPE CanonicalDimension = 3;
    static constexpr UNSIGNED_INTEGER_TYPE Dimension = 3;
    static constexpr UNSIGNED_INTEGER_TYPE NumberOfNodesAtCompileTime = 20;
    static constexpr UNSIGNED_INTEGER_TYPE NumberOfGaussNodesAtCompileTime = 4;

    using BoundaryElementType = Quad<3, Quadratic>;
    static constexpr UNSIGNED_INTEGER_TYPE NumberOfBoundaryElementsAtCompileTime = 6;
};

/**
 * Quadratic Hexahedron (20 nodes)
 *
 *         v
 *         ^
 *         |
 *  3-----10---2
 *  |\     |   |\
 *  | 19   |   | 18
 * 11  \   |   9  \
 *  |   7----14+---6
 *  |   |  +-- |-- | ---> u
 *  0---+-8-\--1   |
 *   \  15   \  \  13
 *   16 |     \  17|
 *     \|      w  \|
 *      4----12----5   u
 *
 */
template<>
struct Hexahedron <Quadratic> : public BaseHexahedron<Hexahedron <Quadratic>> {
    // Types
    using Base = BaseHexahedron<Hexahedron <Quadratic>>;
    using LocalCoordinates = typename Base::LocalCoordinates;
    using WorldCoordinates = typename Base::WorldCoordinates;

    using GaussNode = typename Base::GaussNode;

    template <UNSIGNED_INTEGER_TYPE Dim>
    using Vector = typename Base::template Vector<Dim>;

    template <UNSIGNED_INTEGER_TYPE Rows, UNSIGNED_INTEGER_TYPE Cols>
    using Matrix = typename Base::template Matrix<Rows, Cols>;

    // Constants
    static constexpr auto CanonicalDimension = Base::CanonicalDimension;
    static constexpr auto Dimension = Base::Dimension;
    static constexpr auto NumberOfNodesAtCompileTime = Base::NumberOfNodesAtCompileTime;
    static constexpr auto NumberOfGaussNodesAtCompileTime = Base::NumberOfGaussNodesAtCompileTime;

    // Constructors
    using Base::Base;
    Hexahedron() : Base() {
        this->p_nodes.row(0) = WorldCoordinates(-1, -1, -1); // Node 0
        this->p_nodes.row(1) = WorldCoordinates(+1, -1, -1); // Node 1
        this->p_nodes.row(2) = WorldCoordinates(+1, +1, -1); // Node 2
        this->p_nodes.row(3) = WorldCoordinates(-1, +1, -1); // Node 3
        this->p_nodes.row(4) = WorldCoordinates(-1, -1, +1); // Node 4
        this->p_nodes.row(5) = WorldCoordinates(+1, -1, +1); // Node 5
        this->p_nodes.row(6) = WorldCoordinates(+1, +1, +1); // Node 6
        this->p_nodes.row(7) = WorldCoordinates(-1, +1, +1); // Node 7
        this->p_nodes.row(8)  = WorldCoordinates( 0, -1, -1); // Node 8
        this->p_nodes.row(9)  = WorldCoordinates(+1,  0, -1); // Node 11
        this->p_nodes.row(10) = WorldCoordinates( 0, +1, -1); // Node 13
        this->p_nodes.row(11) = WorldCoordinates(-1,  0, -1); // Node 9
        this->p_nodes.row(12) = WorldCoordinates( 0, -1, +1); // Node 16
        this->p_nodes.row(13) = WorldCoordinates(+1,  0, +1); // Node 18
        this->p_nodes.row(14) = WorldCoordinates( 0, +1, +1); // Node 19
        this->p_nodes.row(15) = WorldCoordinates(-1,  0, +1); // Node 17
        this->p_nodes.row(16) = WorldCoordinates(-1, -1,  0); // Node 10
        this->p_nodes.row(17) = WorldCoordinates(+1, -1,  0); // Node 12
        this->p_nodes.row(18) = WorldCoordinates(+1, +1,  0); // Node 14
        this->p_nodes.row(19) = WorldCoordinates(-1, +1,  0); // Node 15

    };

    // Construct a quadratic Hexahedron from a linear one
    Hexahedron(const Hexahedron<Linear> & linear_Hexahedron) : Base() {
        this->p_nodes.row(0) = linear_Hexahedron.node(0); // Node 0
        this->p_nodes.row(1) = linear_Hexahedron.node(1); // Node 1
        this->p_nodes.row(2) = linear_Hexahedron.node(2); // Node 2
        this->p_nodes.row(3) = linear_Hexahedron.node(3); // Node 3
        this->p_nodes.row(4) = linear_Hexahedron.node(4); // Node 4
        this->p_nodes.row(5) = linear_Hexahedron.node(5); // Node 5
        this->p_nodes.row(6) = linear_Hexahedron.node(6); // Node 6
        this->p_nodes.row(7) = linear_Hexahedron.node(7); // Node 7
        this->p_nodes.row(8)  = linear_Hexahedron.T(LocalCoordinates( 0, -1, -1)); // Node 8
        this->p_nodes.row(9)  = linear_Hexahedron.T(LocalCoordinates(+1,  0, -1)); // Node 9
        this->p_nodes.row(10) = linear_Hexahedron.T(LocalCoordinates( 0, +1, -1)); // Node 10
        this->p_nodes.row(11) = linear_Hexahedron.T(LocalCoordinates(-1,  0, -1)); // Node 11
        this->p_nodes.row(12) = linear_Hexahedron.T(LocalCoordinates( 0, -1, +1)); // Node 12
        this->p_nodes.row(13) = linear_Hexahedron.T(LocalCoordinates(+1,  0, +1)); // Node 13
        this->p_nodes.row(14) = linear_Hexahedron.T(LocalCoordinates( 0, +1, +1)); // Node 14
        this->p_nodes.row(15) = linear_Hexahedron.T(LocalCoordinates(-1,  0, +1)); // Node 15
        this->p_nodes.row(16) = linear_Hexahedron.T(LocalCoordinates(-1, -1,  0)); // Node 16
        this->p_nodes.row(17) = linear_Hexahedron.T(LocalCoordinates(+1, -1,  0)); // Node 17
        this->p_nodes.row(18) = linear_Hexahedron.T(LocalCoordinates(+1, +1,  0)); // Node 18
        this->p_nodes.row(19) = linear_Hexahedron.T(LocalCoordinates(-1, +1,  0)); // Node 19
    }


private:
    // Implementations
    friend struct Element<Hexahedron <Quadratic>>;
    friend struct BaseHexahedron<Hexahedron <Quadratic>>;
    inline auto get_L(const LocalCoordinates & xi) const {
        const auto  & u = xi[0];
        const auto  & v = xi[1];
        const auto  & w = xi[2];

        Vector<NumberOfNodesAtCompileTime> m;
        m[ 8] = 1/4.*(1 - u*u)*(1 - v)*(1 - w);
        m[ 9] = 1/4.*(1 - v*v)*(1 + u)*(1 - w);
        m[10] = 1/4.*(1 - u*u)*(1 + v)*(1 - w);
        m[11] = 1/4.*(1 - v*v)*(1 - u)*(1 - w);
        m[12] = 1/4.*(1 - u*u)*(1 - v)*(1 + w);
        m[13] = 1/4.*(1 - v*v)*(1 + u)*(1 + w);
        m[14] = 1/4.*(1 - u*u)*(1 + v)*(1 + w);
        m[15] = 1/4.*(1 - v*v)*(1 - u)*(1 + w);
        m[16] = 1/4.*(1 - w*w)*(1 - u)*(1 - v);
        m[17] = 1/4.*(1 - w*w)*(1 + u)*(1 - v);
        m[18] = 1/4.*(1 - w*w)*(1 + u)*(1 + v);
        m[19] = 1/4.*(1 - w*w)*(1 - u)*(1 + v);

        m[0] = 1/8.*(1 - u)*(1 - v)*(1 - w) - 1/2.*(m[ 8] + m[11] + m[16]);
        m[1] = 1/8.*(1 + u)*(1 - v)*(1 - w) - 1/2.*(m[ 8] + m[ 9] + m[17]);
        m[2] = 1/8.*(1 + u)*(1 + v)*(1 - w) - 1/2.*(m[ 9] + m[10] + m[18]);
        m[3] = 1/8.*(1 - u)*(1 + v)*(1 - w) - 1/2.*(m[10] + m[11] + m[19]);
        m[4] = 1/8.*(1 - u)*(1 - v)*(1 + w) - 1/2.*(m[12] + m[15] + m[16]);
        m[5] = 1/8.*(1 + u)*(1 - v)*(1 + w) - 1/2.*(m[12] + m[13] + m[17]);
        m[6] = 1/8.*(1 + u)*(1 + v)*(1 + w) - 1/2.*(m[13] + m[14] + m[18]);
        m[7] = 1/8.*(1 - u)*(1 + v)*(1 + w) - 1/2.*(m[14] + m[15] + m[19]);
        return m;
    };

    inline auto get_dL(const LocalCoordinates & xi) const {
        using ShapeDerivative = Vector<3>;
        const auto & u = xi[0];
        const auto & v = xi[1];
        const auto  & w = xi[2];

        Matrix<NumberOfNodesAtCompileTime, CanonicalDimension> m;

        //                                          dL/du                       dL/dv                      dL/dw
        m.row( 8) = ShapeDerivative({-1/2.*u*(1 - v)*(1 - w), -1/4.*(1 - u*u)*(1 - w), -1/4.*(1 - u*u)*(1 - v)});
        m.row( 9) = ShapeDerivative({ 1/4.*(1 - v*v)*(1 - w), -1/2.*v*(1 + u)*(1 - w), -1/4.*(1 - v*v)*(1 + u)});
        m.row(10) = ShapeDerivative({-1/2.*u*(1 + v)*(1 - w),  1/4.*(1 - u*u)*(1 - w), -1/4.*(1 - u*u)*(1 + v)});
        m.row(11) = ShapeDerivative({-1/4.*(1 - v*v)*(1 - w), -1/2.*v*(1 - u)*(1 - w), -1/4.*(1 - v*v)*(1 - u)});
        m.row(12) = ShapeDerivative({-1/2.*u*(1 - v)*(1 + w), -1/4.*(1 - u*u)*(1 + w),  1/4.*(1 - u*u)*(1 - v)});
        m.row(13) = ShapeDerivative({ 1/4.*(1 - v*v)*(1 + w), -1/2.*v*(1 + u)*(1 + w),  1/4.*(1 - v*v)*(1 + u)});
        m.row(14) = ShapeDerivative({-1/2.*u*(1 + v)*(1 + w),  1/4.*(1 - u*u)*(1 + w),  1/4.*(1 - u*u)*(1 + v)});
        m.row(15) = ShapeDerivative({-1/4.*(1 - v*v)*(1 + w), -1/2.*v*(1 - u)*(1 + w),  1/4.*(1 - v*v)*(1 - u)});
        m.row(16) = ShapeDerivative({-1/4.*(1 - w*w)*(1 - v), -1/4.*(1 - w*w)*(1 - u), -1/2.*w*(1 - u)*(1 - v)});
        m.row(17) = ShapeDerivative({ 1/4.*(1 - w*w)*(1 - v), -1/4.*(1 - w*w)*(1 + u), -1/2.*w*(1 + u)*(1 - v)});
        m.row(18) = ShapeDerivative({ 1/4.*(1 - w*w)*(1 + v),  1/4.*(1 - w*w)*(1 + u), -1/2.*w*(1 + u)*(1 + v)});
        m.row(19) = ShapeDerivative({-1/4.*(1 - w*w)*(1 + v),  1/4.*(1 - w*w)*(1 - u), -1/2.*w*(1 - u)*(1 + v)});

        const auto du = m.col(0); // dL/du
        const auto dv = m.col(1); // dL/dv
        const auto dw = m.col(2); // dL/dw
        //                                                        dL/du                                                    dL/dv                                                           dL/dw
        m.row(0) = ShapeDerivative({-1/8.*(1 - v)*(1 - w) - 1/2.*(du[8 ] + du[11] + du[16]), -1/8.*(1 - u)*(1 - w) - 1/2.*(dv[ 8] + dv[11] + dv[16]), -1/8.*(1 - u)*(1 - v) - 1/2.*(dw[ 8] + dw[11] + dw[16])});
        m.row(1) = ShapeDerivative({ 1/8.*(1 - v)*(1 - w) - 1/2.*(du[8 ] + du[ 9] + du[17]), -1/8.*(1 + u)*(1 - w) - 1/2.*(dv[ 8] + dv[ 9] + dv[17]), -1/8.*(1 + u)*(1 - v) - 1/2.*(dw[ 8] + dw[ 9] + dw[17])});
        m.row(2) = ShapeDerivative({ 1/8.*(1 + v)*(1 - w) - 1/2.*(du[9 ] + du[10] + du[18]),  1/8.*(1 + u)*(1 - w) - 1/2.*(dv[ 9] + dv[10] + dv[18]), -1/8.*(1 + u)*(1 + v) - 1/2.*(dw[ 9] + dw[10] + dw[18])});
        m.row(3) = ShapeDerivative({-1/8.*(1 + v)*(1 - w) - 1/2.*(du[10] + du[11] + du[19]),  1/8.*(1 - u)*(1 - w) - 1/2.*(dv[10] + dv[11] + dv[19]), -1/8.*(1 - u)*(1 + v) - 1/2.*(dw[10] + dw[11] + dw[19])});
        m.row(4) = ShapeDerivative({-1/8.*(1 - v)*(1 + w) - 1/2.*(du[12] + du[15] + du[16]), -1/8.*(1 - u)*(1 + w) - 1/2.*(dv[12] + dv[15] + dv[16]),  1/8.*(1 - u)*(1 - v) - 1/2.*(dw[12] + dw[15] + dw[16])});
        m.row(5) = ShapeDerivative({ 1/8.*(1 - v)*(1 + w) - 1/2.*(du[12] + du[13] + du[17]), -1/8.*(1 + u)*(1 + w) - 1/2.*(dv[12] + dv[13] + dv[17]),  1/8.*(1 + u)*(1 - v) - 1/2.*(dw[12] + dw[13] + dw[17])});
        m.row(6) = ShapeDerivative({ 1/8.*(1 + v)*(1 + w) - 1/2.*(du[13] + du[14] + du[18]),  1/8.*(1 + u)*(1 + w) - 1/2.*(dv[13] + dv[14] + dv[18]),  1/8.*(1 + u)*(1 + v) - 1/2.*(dw[13] + dw[14] + dw[18])});
        m.row(7) = ShapeDerivative({-1/8.*(1 + v)*(1 + w) - 1/2.*(du[14] + du[15] + du[19]),  1/8.*(1 - u)*(1 + w) - 1/2.*(dv[14] + dv[15] + dv[19]),  1/8.*(1 - u)*(1 + v) - 1/2.*(dw[14] + dw[15] + dw[19])});


        return m;
    };

    inline auto get_gauss_nodes() const -> const auto & {
        using Weight = FLOATING_POINT_TYPE;
        static const std::vector<GaussNode> gauss_nodes {
            GaussNode {LocalCoordinates(-1/sqrt(3), -1/sqrt(3), -1/sqrt(3)), Weight(1)}, // Node 0
            GaussNode {LocalCoordinates(+1/sqrt(3), -1/sqrt(3), -1/sqrt(3)), Weight(1)}, // Node 1
            GaussNode {LocalCoordinates(-1/sqrt(3), +1/sqrt(3), -1/sqrt(3)), Weight(1)}, // Node 2
            GaussNode {LocalCoordinates(+1/sqrt(3), +1/sqrt(3), -1/sqrt(3)), Weight(1)}, // Node 3
            GaussNode {LocalCoordinates(-1/sqrt(3), -1/sqrt(3), +1/sqrt(3)), Weight(1)}, // Node 4
            GaussNode {LocalCoordinates(+1/sqrt(3), -1/sqrt(3), +1/sqrt(3)), Weight(1)}, // Node 5
            GaussNode {LocalCoordinates(-1/sqrt(3), +1/sqrt(3), +1/sqrt(3)), Weight(1)}, // Node 6
            GaussNode {LocalCoordinates(+1/sqrt(3), +1/sqrt(3), +1/sqrt(3)), Weight(1)}  // Node 7
        };
        return gauss_nodes;
    }

    inline auto get_boundary_elements_nodes() const -> const auto & {
        static const std::array<std::array<UNSIGNED_INTEGER_TYPE, 8>, 6> indices = {{
            {0, 3, 2, 1,  9, 13, 11,  8}, // Face 0
            {0, 4, 7, 3, 10, 17, 15,  9}, // Face 1
            {1, 2, 6, 5, 11, 14, 18, 12}, // Face 2
            {0, 1, 5, 4,  8, 12, 16, 10}, // Face 3
            {2, 3, 7, 6, 13, 15, 19, 14}, // Face 4
            {4, 5, 6, 7, 16, 18, 19, 17}  // Face 5
        }};
        return indices;
    }
};

}
