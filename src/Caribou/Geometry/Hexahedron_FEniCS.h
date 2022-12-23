#pragma once

#include <Caribou/config.h>
#include <Caribou/Geometry/Base/BaseHexahedron.h>
#include <Caribou/Geometry/Quad.h>
#include <Eigen/Core>

namespace caribou::geometry {

struct Hexahedron_FEniCS;

template<>
struct traits<Hexahedron_FEniCS> {
    static constexpr UNSIGNED_INTEGER_TYPE CanonicalDimension = 3;
    static constexpr UNSIGNED_INTEGER_TYPE Dimension = 3;
    static constexpr INTEGER_TYPE NumberOfNodesAtCompileTime = 8;
    static constexpr INTEGER_TYPE NumberOfGaussNodesAtCompileTime = 8;

    using BoundaryElementType = Quad<3>;
    static constexpr INTEGER_TYPE NumberOfBoundaryElementsAtCompileTime = 6;
};

/**
 * Linear Hexahedron_FEniCS
 *
 * \verbatim
 *        v
 * 6----------7
 * |\     ^   |\
 * | \    |   | \
 * |  \   |   |  \
 * |   4------+---5
 * |   |  +-- |-- | -> u
 * 2---+---\--3   |
 *  \  |    \  \  |
 *   \ |     \  \ |
 *    \|      w  \|
 *     0----------1
 * \endverbatim
 *
 */
struct Hexahedron_FEniCS : public BaseHexahedron<Hexahedron_FEniCS> {
    // Types
    using Base = BaseHexahedron<Hexahedron_FEniCS>;
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
    //    u,  v,  w
        {-1, -1, +1}, // Node 0
        {+1, -1, +1}, // Node 1
        {-1, -1, -1}, // Node 2
        {+1, -1, -1}, // Node 3
        {-1, +1, +1}, // Node 4
        {+1, +1, +1}, // Node 5
        {-1, +1, -1}, // Node 6
        {+1, +1, -1}, // Node 7
    };

    // Constructors
    using Base::Base;
    Hexahedron_FEniCS() : Base() {
        this->p_nodes.row(0) = WorldCoordinates(canonical_nodes[0][0], canonical_nodes[0][1], canonical_nodes[0][2]);
        this->p_nodes.row(1) = WorldCoordinates(canonical_nodes[1][0], canonical_nodes[1][1], canonical_nodes[1][2]);
        this->p_nodes.row(2) = WorldCoordinates(canonical_nodes[2][0], canonical_nodes[2][1], canonical_nodes[2][2]);
        this->p_nodes.row(3) = WorldCoordinates(canonical_nodes[3][0], canonical_nodes[3][1], canonical_nodes[3][2]);
        this->p_nodes.row(4) = WorldCoordinates(canonical_nodes[4][0], canonical_nodes[4][1], canonical_nodes[4][2]);
        this->p_nodes.row(5) = WorldCoordinates(canonical_nodes[5][0], canonical_nodes[5][1], canonical_nodes[5][2]);
        this->p_nodes.row(6) = WorldCoordinates(canonical_nodes[6][0], canonical_nodes[6][1], canonical_nodes[6][2]);
        this->p_nodes.row(7) = WorldCoordinates(canonical_nodes[7][0], canonical_nodes[7][1], canonical_nodes[7][2]);
    };

    // Public methods
    /**
     * Get the list of node indices of the edges.
     * \sa Element::boundary_elements_node_indices
     */
    inline auto edges() const {
        static const std::array<std::array<UNSIGNED_INTEGER_TYPE, 2>, 12> indices = {{
            {2, 3}, // Edge 0
            {3, 7}, // Edge 1
            {6, 7}, // Edge 2
            {6, 2}, // Edge 3
            {0, 1}, // Edge 4
            {1, 5}, // Edge 5
            {5, 4}, // Edge 6
            {4, 0}, // Edge 7
            {2, 0}, // Edge 8
            {3, 1}, // Edge 9
            {7, 5}, // Edge 10
            {6, 4}  // Edge 11
        }};
        return indices;
    }


private:
    // Implementations
    friend struct Element<Hexahedron_FEniCS>;
    friend struct BaseHexahedron<Hexahedron_FEniCS>;
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
            {2, 6, 7, 3}, // Face 0
            {2, 0, 4, 6}, // Face 1
            {3, 7, 5, 1}, // Face 2
            {2, 3, 1, 0}, // Face 3
            {7, 6, 4, 5}, // Face 4
            {0, 1, 5, 4}  // Face 5
        }};
        return indices;
    }
};
}