#pragma once

#include <Caribou/config.h>
#include <Caribou/Geometry/Base/BaseQuad.h>
#include <Caribou/Geometry/Segment.h>
#include <Eigen/Core>
#include <array>

namespace caribou::geometry {

template<UNSIGNED_INTEGER_TYPE Dimension>
struct Quad;

template<UNSIGNED_INTEGER_TYPE _Dimension>
struct traits<Quad <_Dimension>> {
    static constexpr UNSIGNED_INTEGER_TYPE CanonicalDimension = 2;
    static constexpr UNSIGNED_INTEGER_TYPE Dimension = _Dimension;
    static constexpr INTEGER_TYPE NumberOfNodesAtCompileTime = 4;
    static constexpr INTEGER_TYPE NumberOfGaussNodesAtCompileTime = 4;

    using BoundaryElementType = Segment<_Dimension>;
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
struct Quad: public BaseQuad<Quad <_Dimension>> {
    // Types
    using Base = BaseQuad<Quad <_Dimension>>;
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
    friend struct Element<Quad <_Dimension>>;
    friend struct BaseQuad<Quad <_Dimension>>;
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

}
