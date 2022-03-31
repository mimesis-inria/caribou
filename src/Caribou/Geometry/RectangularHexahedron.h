#pragma once

#include <Caribou/config.h>
#include <Caribou/Geometry/Base/BaseRectangularHexahedron.h>
#include <Caribou/Geometry/Hexahedron.h>
#include <Caribou/Geometry/RectangularQuad.h>
#include <Eigen/Core>

namespace caribou::geometry {

struct RectangularHexahedron;

template<>
struct traits<RectangularHexahedron> {
    static constexpr UNSIGNED_INTEGER_TYPE CanonicalDimension = 3;
    static constexpr UNSIGNED_INTEGER_TYPE Dimension = 3;
    static constexpr INTEGER_TYPE NumberOfNodesAtCompileTime = 8;
    static constexpr INTEGER_TYPE NumberOfGaussNodesAtCompileTime = 8;

    using BoundaryElementType = RectangularQuad<3>;
    static constexpr INTEGER_TYPE NumberOfBoundaryElementsAtCompileTime = 6;
};

/**
 * Linear Rectangular Hexahedron
 *
 * \verbatim
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
 * \endverbatim
 *
 */
struct RectangularHexahedron: public BaseRectangularHexahedron<RectangularHexahedron> {
    // Types
    using Base = BaseRectangularHexahedron<RectangularHexahedron>;
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

    static constexpr auto canonical_nodes = Hexahedron::canonical_nodes;

    // Constructors
    using Base::Base;

    // Constructor from a regular Hexahedron
    explicit RectangularHexahedron(const Hexahedron & hexa) {
        this->p_center = hexa.center();
        this->p_H = Size((hexa.node(1)-hexa.node(0)).norm(), (hexa.node(3)-hexa.node(0)).norm(), (hexa.node(4)-hexa.node(0)).norm());
        this->p_R = hexa.frame({0, 0, 0});
    };

    /** Constructor from an Eigen matrix containing the positions of the hexa's nodes */
    template<typename EigenType, REQUIRES(EigenType::RowsAtCompileTime == NumberOfNodesAtCompileTime)>
    explicit RectangularHexahedron(Eigen::EigenBase<EigenType> & nodes)
        : RectangularHexahedron(Hexahedron(nodes)) {}

    /** Constructor from an Eigen matrix containing the positions of the hexa's nodes */
    template<typename EigenType, REQUIRES(EigenType::RowsAtCompileTime == NumberOfNodesAtCompileTime)>
    explicit RectangularHexahedron(const Eigen::EigenBase<EigenType> & nodes)
        : RectangularHexahedron(Hexahedron(nodes)) {}

    // Public methods
    /**
     * Get the list of node indices of the edges.
     * \sa Element::boundary_elements_node_indices
     */
    inline auto edges() const {
        static const std::array<std::array<UNSIGNED_INTEGER_TYPE, 2>, 12> indices = {{
            {0, 1}, // Edge 0
            {1, 2}, // Edge 1
            {2, 3}, // Edge 2
            {3, 0}, // Edge 3
            {4, 5}, // Edge 4
            {5, 6}, // Edge 5
            {6, 7}, // Edge 6
            {7, 4}, // Edge 7
            {0, 4}, // Edge 8
            {1, 5}, // Edge 9
            {2, 6}, // Edge 10
            {3, 7}  // Edge 11
        }};
        return indices;
    }

private:
    // Implementations
    friend struct Element<RectangularHexahedron>;
    friend struct BaseRectangularHexahedron<RectangularHexahedron>;
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

    inline auto get_node(const UNSIGNED_INTEGER_TYPE & index) const -> WorldCoordinates {
        caribou_assert(index < NumberOfNodesAtCompileTime and "Trying to get a node from an invalid node index.");
        return this->world_coordinates(LocalCoordinates(canonical_nodes[index][0], canonical_nodes[index][1], canonical_nodes[index][2]));
    }

    inline auto get_nodes() const {
        Matrix<NumberOfNodesAtCompileTime, Dimension> nodes;
        for (UNSIGNED_INTEGER_TYPE node_id = 0; node_id < NumberOfNodesAtCompileTime; ++node_id) {
            nodes.row(node_id) = get_node(node_id);
        }
        return nodes;
    }

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

}
