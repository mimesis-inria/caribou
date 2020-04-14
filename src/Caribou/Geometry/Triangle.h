#pragma once

#include <Caribou/config.h>
#include <Caribou/Geometry/BaseTriangle.h>
#include <Caribou/Geometry/Segment.h>
#include <Eigen/Core>
#include <array>

namespace caribou::geometry {

template<UNSIGNED_INTEGER_TYPE _Dimension, UNSIGNED_INTEGER_TYPE _Order = Linear>
struct Triangle;

template<UNSIGNED_INTEGER_TYPE _Dimension>
struct traits<Triangle <_Dimension, Linear>> {
    static constexpr UNSIGNED_INTEGER_TYPE CanonicalDimension = 2;
    static constexpr UNSIGNED_INTEGER_TYPE Dimension = _Dimension;
    static constexpr INTEGER_TYPE NumberOfNodesAtCompileTime = 3;
    static constexpr INTEGER_TYPE NumberOfGaussNodesAtCompileTime = 1;

    using BoundaryElementType = Segment<_Dimension, Linear>;
    static constexpr INTEGER_TYPE NumberOfBoundaryElementsAtCompileTime = 3;
};

/**
 * Linear Triangle
 *
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
 *
 * @tparam _Dimension The world coordinates dimension
 */
template<UNSIGNED_INTEGER_TYPE _Dimension>
struct Triangle <_Dimension, Linear> : public BaseTriangle<Triangle <_Dimension, Linear>> {
    // Types
    using Base = BaseTriangle<Triangle <_Dimension, Linear>>;
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
    friend struct Element<Triangle <_Dimension, Linear>>;
    friend struct BaseTriangle<Triangle <_Dimension, Linear>>;
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


template<UNSIGNED_INTEGER_TYPE _Dimension>
struct traits<Triangle <_Dimension, Quadratic>> {
    static constexpr UNSIGNED_INTEGER_TYPE CanonicalDimension = 2;
    static constexpr UNSIGNED_INTEGER_TYPE Dimension = _Dimension;
    static constexpr INTEGER_TYPE NumberOfNodesAtCompileTime = 6;
    static constexpr INTEGER_TYPE NumberOfGaussNodesAtCompileTime = 3;

    using BoundaryElementType = Segment<_Dimension, Quadratic>;
    static constexpr INTEGER_TYPE NumberOfBoundaryElementsAtCompileTime = 3;
};

/**
 * Quadratic Triangle
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
 * @tparam _Dimension The world coordinates dimension
 */
template<UNSIGNED_INTEGER_TYPE _Dimension>
struct Triangle <_Dimension, Quadratic> : public BaseTriangle<Triangle <_Dimension, Quadratic>> {
    // Types
    using Base = BaseTriangle<Triangle <_Dimension, Quadratic>>;
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
            this->p_nodes.row(0) = WorldCoordinates(0.0, 0.0); // Node 0
            this->p_nodes.row(1) = WorldCoordinates(1.0, 0.0); // Node 1
            this->p_nodes.row(2) = WorldCoordinates(0.0, 1.0); // Node 2
            this->p_nodes.row(3) = WorldCoordinates(0.5, 0.0); // Node 3
            this->p_nodes.row(4) = WorldCoordinates(0.5, 0.5); // Node 4
            this->p_nodes.row(5) = WorldCoordinates(0.0, 0.5); // Node 5
        } else {
            this->p_nodes.row(0) = WorldCoordinates(0.0, 0.0, 0.0); // Node 0
            this->p_nodes.row(1) = WorldCoordinates(1.0, 0.0, 0.0); // Node 1
            this->p_nodes.row(2) = WorldCoordinates(0.0, 1.0, 0.0); // Node 2
            this->p_nodes.row(3) = WorldCoordinates(0.5, 0.0, 0.0); // Node 3
            this->p_nodes.row(4) = WorldCoordinates(0.5, 0.5, 0.0); // Node 4
            this->p_nodes.row(5) = WorldCoordinates(0.0, 0.5, 0.0); // Node 5
        }
    };

    // Construct a quadratic triangle from a linear one
    explicit Triangle(const Triangle<Dimension, Linear> & linear_triangle) : Base() {
        this->p_nodes.row(0) = linear_triangle.node(0); // Node 0
        this->p_nodes.row(1) = linear_triangle.node(1); // Node 1
        this->p_nodes.row(2) = linear_triangle.node(2); // Node 2
        this->p_nodes.row(3) = linear_triangle.world_coordinates(LocalCoordinates(0.5, 0.0)); // Node 3
        this->p_nodes.row(4) = linear_triangle.world_coordinates(LocalCoordinates(0.5, 0.5)); // Node 4
        this->p_nodes.row(5) = linear_triangle.world_coordinates(LocalCoordinates(0.0, 0.5)); // Node 5
    }

    /** Construct a quadratic triangle from three nodes */
    Triangle(WorldCoordinates & p0, WorldCoordinates & p1, WorldCoordinates & p2)
    : Triangle(Triangle<Dimension, Linear>(p0, p1, p2)) {}

    /** Construct a quadratic triangle from three nodes */
    Triangle(const WorldCoordinates & p0, const WorldCoordinates & p1, const WorldCoordinates & p2)
        : Triangle(Triangle<Dimension, Linear>(p0, p1, p2)) {}

private:
    // Implementations
    friend struct Element<Triangle <_Dimension, Quadratic>>;
    friend struct BaseTriangle<Triangle <_Dimension, Quadratic>>;
    inline auto get_L(const LocalCoordinates & xi) const {
        const auto  & u = xi[0];
        const auto  & v = xi[1];
        const auto l = 1 - u - v;

        Vector<NumberOfNodesAtCompileTime> m;
        m << l * (2*l - 1),   // Node 0
             u * (2*u - 1),   // Node 1
             v * (2*v - 1),   // Node 2
             4 * u * l,       // Node 3
             4 * u * v,       // Node 4
             4 * v * l;       // Node 5

        return m;
    };

    inline auto get_dL(const LocalCoordinates & xi) const {
        const auto & u = xi[0];
        const auto & v = xi[1];
        const auto l = 1 - u - v;

        Matrix<NumberOfNodesAtCompileTime, CanonicalDimension> m;
        //       dL/du              dL/dv
        m <<   1 - 4 * l  ,       1 - 4 * l  ,   // Node 0
               4 * u - 1  ,          0       ,   // Node 1
                    0     ,       4 * v - 1  ,   // Node 2
               4 * (l - u),      -4 * u      ,   // Node 3
                  4 * v   ,       4 * u      ,   // Node 4
                - 4 * v   ,       4 * (l - v);   // Node 5
        return m;
    };

    inline auto get_gauss_nodes() const -> const auto & {
        using Weight = FLOATING_POINT_TYPE;
        static const std::vector<GaussNode> gauss_nodes {
            GaussNode {LocalCoordinates(2/3., 1/6.), Weight(1/6.)}, // Node 0
            GaussNode {LocalCoordinates(1/6., 2/3.), Weight(1/6.)}, // Node 1
            GaussNode {LocalCoordinates(1/6., 1/6.), Weight(1/6.)}  // Node 2
        };
        return gauss_nodes;
    }

    inline auto get_boundary_elements_nodes() const -> const auto & {
        static const std::array<std::array<UNSIGNED_INTEGER_TYPE, 3>, 3> indices = {{
            {0, 1, 3}, // Edge 0
            {1, 2, 4}, // Edge 1
            {2, 0, 5}  // Edge 2
        }};
        return indices;
    }
};

}
