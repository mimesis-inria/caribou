#pragma once

#include <Caribou/config.h>
#include <Caribou/Geometry/BaseQuad.h>
#include <Caribou/Geometry/Segment.h>
#include <Eigen/Core>
#include <array>

namespace caribou::geometry {

template<UNSIGNED_INTEGER_TYPE Dimension, UNSIGNED_INTEGER_TYPE Order = Linear>
struct Quad;

template<UNSIGNED_INTEGER_TYPE _Dimension>
struct traits<Quad <_Dimension, Linear>> {
    static constexpr UNSIGNED_INTEGER_TYPE CanonicalDimension = 2;
    static constexpr UNSIGNED_INTEGER_TYPE Dimension = _Dimension;
    static constexpr INTEGER_TYPE NumberOfNodesAtCompileTime = 4;
    static constexpr INTEGER_TYPE NumberOfGaussNodesAtCompileTime = 4;

    using BoundaryElementType = Segment<_Dimension, Linear>;
    static constexpr INTEGER_TYPE NumberOfBoundaryElementsAtCompileTime = 4;
};

/**
 * Linear Quad
 *
 * \verbatim
 *        v
 *        ^
 *        |
 *  3-----------2
 *  |     |     |
 *  |     |     |
 *  |     +---- | --> u
 *  |           |
 *  |           |
 *  0-----------1
 *  \endverbatim
 *
 * @tparam _Dimension The world coordinates dimension
 */
template<UNSIGNED_INTEGER_TYPE _Dimension>
struct Quad <_Dimension, Linear> : public BaseQuad<Quad <_Dimension, Linear>> {
    // Types
    using Base = BaseQuad<Quad <_Dimension, Linear>>;
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

    static constexpr FLOATING_POINT_TYPE canonical_nodes [NumberOfNodesAtCompileTime][CanonicalDimension] = {
    //    u,  v
        {-1, -1}, // Node 0
        {+1, -1}, // Node 1
        {+1, +1}, // Node 2
        {-1, +1}  // Node 3
    };

    // Constructors
    using Base::Base;
    Quad() : Base() {
        if constexpr (Dimension == 2) {
            this->p_nodes.row(0) = WorldCoordinates(canonical_nodes[0][0], canonical_nodes[0][1]); // Node 0
            this->p_nodes.row(1) = WorldCoordinates(canonical_nodes[1][0], canonical_nodes[1][1]); // Node 1
            this->p_nodes.row(2) = WorldCoordinates(canonical_nodes[2][0], canonical_nodes[2][1]); // Node 2
            this->p_nodes.row(3) = WorldCoordinates(canonical_nodes[3][0], canonical_nodes[3][1]); // Node 3
        } else {
            this->p_nodes.row(0) = WorldCoordinates(canonical_nodes[0][0], canonical_nodes[0][1], 0); // Node 0
            this->p_nodes.row(1) = WorldCoordinates(canonical_nodes[1][0], canonical_nodes[1][1], 0); // Node 1
            this->p_nodes.row(2) = WorldCoordinates(canonical_nodes[2][0], canonical_nodes[2][1], 0); // Node 2
            this->p_nodes.row(3) = WorldCoordinates(canonical_nodes[3][0], canonical_nodes[3][1], 0); // Node 3
        }
    };


private:
    // Implementations
    friend struct Element<Quad <_Dimension, Linear>>;
    friend struct BaseQuad<Quad <_Dimension, Linear>>;
    inline auto get_L(const LocalCoordinates & xi) const -> Vector<NumberOfNodesAtCompileTime> {
        const auto  & u = xi[0];
        const auto  & v = xi[1];
        return {
            static_cast<FLOATING_POINT_TYPE>(1 / 4. * (1 - u) * (1 - v)),
            static_cast<FLOATING_POINT_TYPE>(1 / 4. * (1 + u) * (1 - v)),
            static_cast<FLOATING_POINT_TYPE>(1 / 4. * (1 + u) * (1 + v)),
            static_cast<FLOATING_POINT_TYPE>(1 / 4. * (1 - u) * (1 + v))
        };
    };

    inline auto get_dL(const LocalCoordinates & xi) const {
        const auto & u = xi[0];
        const auto & v = xi[1];

        Matrix<NumberOfNodesAtCompileTime, CanonicalDimension> m;
        //         dL/du                dL/dv
        m << -1 / 4. * (1 - v),   -1 / 4. * (1 - u),   // Node 0
             +1 / 4. * (1 - v),   -1 / 4. * (1 + u),   // Node 1
             +1 / 4. * (1 + v),   +1 / 4. * (1 + u),   // Node 2
             -1 / 4. * (1 + v),   +1 / 4. * (1 - u);   // Node 3
        return m;
    };

    inline auto get_gauss_nodes() const -> const auto & {
        using Weight = FLOATING_POINT_TYPE;
        static const std::vector<GaussNode> gauss_nodes {
            GaussNode {LocalCoordinates(-1/sqrt(3), -1/sqrt(3)), Weight(1)}, // Node 0
            GaussNode {LocalCoordinates(+1/sqrt(3), -1/sqrt(3)), Weight(1)}, // Node 1
            GaussNode {LocalCoordinates(+1/sqrt(3), +1/sqrt(3)), Weight(1)}, // Node 2
            GaussNode {LocalCoordinates(-1/sqrt(3), +1/sqrt(3)), Weight(1)}  // Node 3
        };
        return gauss_nodes;
    }

    inline auto get_boundary_elements_nodes() const -> const auto & {
        static const std::array<std::array<UNSIGNED_INTEGER_TYPE, 2>, 4> indices = {{
            {0, 1}, // Edge 0
            {1, 2}, // Edge 1
            {2, 3}, // Edge 2
            {3, 0}  // Edge 3
        }};
        return indices;
    }
};


template<UNSIGNED_INTEGER_TYPE _Dimension>
struct traits<Quad <_Dimension, Quadratic>> {
    static constexpr UNSIGNED_INTEGER_TYPE CanonicalDimension = 2;
    static constexpr UNSIGNED_INTEGER_TYPE Dimension = _Dimension;
    static constexpr INTEGER_TYPE NumberOfNodesAtCompileTime = 8;
    static constexpr INTEGER_TYPE NumberOfGaussNodesAtCompileTime = 4;

    using BoundaryElementType = Segment<_Dimension, Quadratic>;
    static constexpr INTEGER_TYPE NumberOfBoundaryElementsAtCompileTime = 4;
};

/**
 * Quadratic Quad
 *
 * \verbatim
 *       v
 *       ^
 *       |
 * 3-----6-----2
 * |     |     |
 * |     |     |
 * 7     +---- 5 --> u
 * |           |
 * |           |
 * 0-----4-----1
 * \endverbatim
 *
 * @tparam _Dimension The world coordinates dimension
 */
template<UNSIGNED_INTEGER_TYPE _Dimension>
struct Quad <_Dimension, Quadratic> : public BaseQuad<Quad <_Dimension, Quadratic>> {
    // Types
    using Base = BaseQuad<Quad <_Dimension, Quadratic>>;
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

    static constexpr FLOATING_POINT_TYPE canonical_nodes [NumberOfNodesAtCompileTime][CanonicalDimension] = {
    //    u,  v
        {-1, -1}, // Node 0
        {+1, -1}, // Node 1
        {+1, +1}, // Node 2
        {-1, +1}, // Node 3
        { 0, -1}, // Node 4
        {+1,  0}, // Node 5
        { 0, +1}, // Node 6
        {-1,  0}  // Node 7
    };

    // Constructors
    using Base::Base;
    Quad() : Base() {
        if constexpr (Dimension == 2) {
            this->p_nodes.row(0) = WorldCoordinates(canonical_nodes[0][0], canonical_nodes[0][1]); // Node 0
            this->p_nodes.row(1) = WorldCoordinates(canonical_nodes[1][0], canonical_nodes[1][1]); // Node 1
            this->p_nodes.row(2) = WorldCoordinates(canonical_nodes[2][0], canonical_nodes[2][1]); // Node 2
            this->p_nodes.row(3) = WorldCoordinates(canonical_nodes[3][0], canonical_nodes[3][1]); // Node 3
            this->p_nodes.row(4) = WorldCoordinates(canonical_nodes[4][0], canonical_nodes[4][1]); // Node 4
            this->p_nodes.row(5) = WorldCoordinates(canonical_nodes[5][0], canonical_nodes[5][1]); // Node 5
            this->p_nodes.row(6) = WorldCoordinates(canonical_nodes[6][0], canonical_nodes[6][1]); // Node 6
            this->p_nodes.row(7) = WorldCoordinates(canonical_nodes[7][0], canonical_nodes[7][1]); // Node 7
        } else {
            this->p_nodes.row(0) = WorldCoordinates(canonical_nodes[0][0], canonical_nodes[0][1], 0); // Node 0
            this->p_nodes.row(1) = WorldCoordinates(canonical_nodes[1][0], canonical_nodes[1][1], 0); // Node 1
            this->p_nodes.row(2) = WorldCoordinates(canonical_nodes[2][0], canonical_nodes[2][1], 0); // Node 2
            this->p_nodes.row(3) = WorldCoordinates(canonical_nodes[3][0], canonical_nodes[3][1], 0); // Node 3
            this->p_nodes.row(4) = WorldCoordinates(canonical_nodes[4][0], canonical_nodes[4][1], 0); // Node 4
            this->p_nodes.row(5) = WorldCoordinates(canonical_nodes[5][0], canonical_nodes[5][1], 0); // Node 5
            this->p_nodes.row(6) = WorldCoordinates(canonical_nodes[6][0], canonical_nodes[6][1], 0); // Node 6
            this->p_nodes.row(7) = WorldCoordinates(canonical_nodes[7][0], canonical_nodes[7][1], 0); // Node 7
        }
    };

    /** Construct a quadratic Quad from a linear one */
    explicit Quad(const Quad<Dimension, Linear> & linear_Quad) : Base() {
        this->p_nodes.row(0) = linear_Quad.node(0); // Node 0
        this->p_nodes.row(1) = linear_Quad.node(1); // Node 1
        this->p_nodes.row(2) = linear_Quad.node(2); // Node 2
        this->p_nodes.row(3) = linear_Quad.node(3); // Node 3
        this->p_nodes.row(4) = linear_Quad.world_coordinates(LocalCoordinates(canonical_nodes[4][0], canonical_nodes[4][1])); // Node 4
        this->p_nodes.row(5) = linear_Quad.world_coordinates(LocalCoordinates(canonical_nodes[5][0], canonical_nodes[5][1])); // Node 5
        this->p_nodes.row(6) = linear_Quad.world_coordinates(LocalCoordinates(canonical_nodes[6][0], canonical_nodes[6][1])); // Node 6
        this->p_nodes.row(7) = linear_Quad.world_coordinates(LocalCoordinates(canonical_nodes[7][0], canonical_nodes[7][1])); // Node 7
    }

    /** Construct a quadratic Quad from four nodes */
    Quad(WorldCoordinates & p0, WorldCoordinates & p1, WorldCoordinates & p2, WorldCoordinates & p3)
    : Quad(Quad<Dimension, Linear>(p0, p1, p2, p3)) {}

    /** Construct a quadratic Quad from four nodes */
    Quad(const WorldCoordinates & p0, const WorldCoordinates & p1, const WorldCoordinates & p2, const WorldCoordinates & p3)
    : Quad(Quad<Dimension, Linear>(p0, p1, p2, p3)) {}


private:
    // Implementations
    friend struct Element<Quad <_Dimension, Quadratic>>;
    friend struct BaseQuad<Quad <_Dimension, Quadratic>>;
    inline auto get_L(const LocalCoordinates & xi) const {
        const auto  & u = xi[0];
        const auto  & v = xi[1];

        Vector<NumberOfNodesAtCompileTime> m;
        m[4] = 1/2.*(1 - u*u)*(1 - v); // Node 4
        m[5] = 1/2.*(1 - v*v)*(1 + u); // Node 5
        m[6] = 1/2.*(1 - u*u)*(1 + v); // Node 6
        m[7] = 1/2.*(1 - v*v)*(1 - u); // Node 7

        m[0] = 1/4.*(1 - u)*(1 - v) - 1/2.*(m[4] + m[7]); // Node 0
        m[1] = 1/4.*(1 + u)*(1 - v) - 1/2.*(m[4] + m[5]); // Node 1
        m[2] = 1/4.*(1 + u)*(1 + v) - 1/2.*(m[5] + m[6]); // Node 2
        m[3] = 1/4.*(1 - u)*(1 + v) - 1/2.*(m[6] + m[7]); // Node 3

        return m;
    };

    inline auto get_dL(const LocalCoordinates & xi) const {
        const auto & u = xi[0];
        const auto & v = xi[1];

        Matrix<NumberOfNodesAtCompileTime, CanonicalDimension> m;
        //                          dL/du                                                        dL/dv
        m <<   -1/4.*(1 - v) - 1/2.*(-u*(1 - v) + -1/2.*(1 - v*v))  ,      -1/4.*(1 - u) - 1/2.*(-1/2.*(1 - u*u) + -v*(1 - u)),   // Node 0
                1/4.*(1 - v) - 1/2.*(-u*(1 - v) +  1/2.*(1 - v*v))  ,      -1/4.*(1 + u) - 1/2.*(-1/2.*(1 - u*u) + -v*(1 + u)),   // Node 1
                1/4.*(1 + v) - 1/2.*(-u*(1 + v) +  1/2.*(1 - v*v))  ,       1/4.*(1 + u) - 1/2.*( 1/2.*(1 - u*u) + -v*(1 + u)),   // Node 2
               -1/4.*(1 + v) - 1/2.*(-u*(1 + v) + -1/2.*(1 - v*v))  ,       1/4.*(1 - u) - 1/2.*( 1/2.*(1 - u*u) + -v*(1 - u)),   // Node 3
                                  -u*(1 - v)                        ,                       -1/2.*(1 - u*u)                   ,   // Node 4
                               1/2.*(1 - v*v)                       ,                         -v*(1 + u)                      ,   // Node 5
                                  -u*(1 + v)                        ,                        1/2.*(1 - u*u)                   ,   // Node 6
                               -1/2.*(1 - v*v)                      ,                         -v*(1 - u);                         // Node 7

        return m;
    };

    inline auto get_gauss_nodes() const -> const auto & {
        using Weight = FLOATING_POINT_TYPE;
        static const std::vector<GaussNode> gauss_nodes {
            GaussNode {LocalCoordinates(-1/sqrt(3), -1/sqrt(3)), Weight(1)}, // Node 0
            GaussNode {LocalCoordinates(+1/sqrt(3), -1/sqrt(3)), Weight(1)}, // Node 1
            GaussNode {LocalCoordinates(+1/sqrt(3), +1/sqrt(3)), Weight(1)}, // Node 2
            GaussNode {LocalCoordinates(-1/sqrt(3), +1/sqrt(3)), Weight(1)}  // Node 3
        };
        return gauss_nodes;
    }

    inline auto get_boundary_elements_nodes() const -> const auto & {
        static const std::array<std::array<UNSIGNED_INTEGER_TYPE, 3>, 4> indices = {{
            {0, 1, 4}, // Edge 0
            {1, 2, 5}, // Edge 1
            {2, 3, 6}, // Edge 2
            {3, 0, 7}  // Edge 3
        }};
        return indices;
    }
};

}
