#pragma once

#include <Caribou/config.h>
#include <Caribou/Geometry/Base/BaseTriangle.h>
#include <Caribou/Geometry/Segment.h>
#include <Eigen/Core>
#include <array>

namespace caribou::geometry {

template<UNSIGNED_INTEGER_TYPE _Dimension>
struct Triangle;

template<UNSIGNED_INTEGER_TYPE _Dimension>
struct traits<Triangle <_Dimension>> {
    static constexpr UNSIGNED_INTEGER_TYPE CanonicalDimension = 2;
    static constexpr UNSIGNED_INTEGER_TYPE Dimension = _Dimension;
    static constexpr INTEGER_TYPE NumberOfNodesAtCompileTime = 3;
    static constexpr INTEGER_TYPE NumberOfGaussNodesAtCompileTime = 1;

    using BoundaryElementType = Segment<_Dimension>;
    static constexpr INTEGER_TYPE NumberOfBoundaryElementsAtCompileTime = 3;
};

/**
 * Linear Triangle
 *
 * \verbatim
 * v
 * ^
 * |
 * 2
 * |`\
 * |  `\
 * |    `\
 * |      `\
 * |        `\
 * 0----------1 --> u
 * \endverbatim
 *
 * @tparam _Dimension The world coordinates dimension
 */
template<UNSIGNED_INTEGER_TYPE _Dimension>
struct Triangle: public BaseTriangle<Triangle <_Dimension>> {
    // Types
    using Base = BaseTriangle<Triangle <_Dimension>>;
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
    Triangle() : Base() {
        if constexpr (Dimension == 2) {
            this->p_nodes.row(0) = WorldCoordinates(0, 0);
            this->p_nodes.row(1) = WorldCoordinates(1, 0);
            this->p_nodes.row(2) = WorldCoordinates(0, 1);
        } else {
            this->p_nodes.row(0) = WorldCoordinates(0, 0, 0);
            this->p_nodes.row(1) = WorldCoordinates(1, 0, 0);
            this->p_nodes.row(2) = WorldCoordinates(0, 1, 0);
        }
    };


private:
    // Implementations
    friend struct Element<Triangle <_Dimension>>;
    friend struct BaseTriangle<Triangle <_Dimension>>;
    inline auto get_L(const LocalCoordinates & xi) const -> Vector<NumberOfNodesAtCompileTime> {
        const auto  & u = xi[0];
        const auto  & v = xi[1];
        return {
            static_cast<FLOATING_POINT_TYPE>(1 - u - v),
            static_cast<FLOATING_POINT_TYPE>(u),
            static_cast<FLOATING_POINT_TYPE>(v)
        };
    };

    inline auto get_dL(const LocalCoordinates & /*xi*/) const {
        Matrix<NumberOfNodesAtCompileTime, CanonicalDimension> m;
        //    dL/du dL/dv
        m <<   -1,   -1,   // Node 0
                1,    0,   // Node 1
                0,    1;   // Node 2
        return m;
    };

    inline auto get_gauss_nodes() const -> const auto & {
        using Weight = FLOATING_POINT_TYPE;
        static const std::vector<GaussNode> gauss_nodes {
            GaussNode {LocalCoordinates(1/3., 1/3.), Weight(1/2.)} // Node 0
        };
        return gauss_nodes;
    }

    inline auto get_boundary_elements_nodes() const -> const auto & {
        static const std::array<std::array<UNSIGNED_INTEGER_TYPE, 2>, 3> indices = {{
            {0, 1}, // Edge 0
            {1, 2}, // Edge 1
            {2, 0}  // Edge 2
        }};
        return indices;
    }
};

}
