#pragma once

#include <Caribou/config.h>
#include <Caribou/Geometry/Base/BaseRectangularQuad.h>
#include <Caribou/Geometry/Quad.h>
#include <Eigen/Core>
#include <array>

namespace caribou::geometry {

template<UNSIGNED_INTEGER_TYPE Dimension>
struct RectangularQuad;


template<UNSIGNED_INTEGER_TYPE _Dimension>
struct traits<RectangularQuad <_Dimension>> {
    static constexpr UNSIGNED_INTEGER_TYPE CanonicalDimension = 2;
    static constexpr UNSIGNED_INTEGER_TYPE Dimension = _Dimension;
    static constexpr INTEGER_TYPE NumberOfNodesAtCompileTime = 4;
    static constexpr INTEGER_TYPE NumberOfGaussNodesAtCompileTime = 4;

    using BoundaryElementType = Segment<_Dimension>;
    static constexpr INTEGER_TYPE NumberOfBoundaryElementsAtCompileTime = 4;
};

/**
 * Linear Rectangular Quad
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
 * \endverbatim
 *
 * @tparam _Dimension The world coordinates dimension
 */
template<UNSIGNED_INTEGER_TYPE _Dimension>
struct RectangularQuad : public BaseRectangularQuad<RectangularQuad <_Dimension>> {
    // Types
    using Base = BaseRectangularQuad<RectangularQuad <_Dimension>>;
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

    static constexpr auto canonical_nodes = Quad<Dimension>::canonical_nodes;

    using Size = typename Base::Size;
    using Rotation = typename Base::Rotation;

    // Constructors
    using Base::Base;

    // Constructor from a regular quad
    explicit RectangularQuad(const Quad<Dimension> & quad) {
        this->p_center = quad.center();
        this->p_H = Size((quad.node(1)-quad.node(0)).norm(), (quad.node(3)-quad.node(0)).norm());
        this->p_R = quad.frame({0, 0});
    };

    /** Constructor from an Eigen matrix containing the positions of the quad's nodes */
    template<typename EigenType, REQUIRES(EigenType::RowsAtCompileTime == NumberOfNodesAtCompileTime)>
    explicit RectangularQuad(Eigen::EigenBase<EigenType> & nodes)
    : RectangularQuad(Quad<Dimension>(nodes)) {}

    /** Constructor from an Eigen matrix containing the positions of the quad's nodes */
    template<typename EigenType, REQUIRES(EigenType::RowsAtCompileTime == NumberOfNodesAtCompileTime)>
    explicit RectangularQuad(const Eigen::EigenBase<EigenType> & nodes)
    : RectangularQuad(Quad<Dimension>(nodes)) {}

private:
    // Implementations
    friend struct Element<RectangularQuad <_Dimension>>;
    friend struct BaseQuad<RectangularQuad <_Dimension>>;
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

    inline auto get_node(const UNSIGNED_INTEGER_TYPE & index) const -> WorldCoordinates {
        caribou_assert(index < NumberOfNodesAtCompileTime and "Trying to get a node from an invalid node index.");
        return this->T(LocalCoordinates(canonical_nodes[index][0], canonical_nodes[index][1]));
    }

    inline auto get_nodes() const {
        Matrix<NumberOfNodesAtCompileTime, Dimension> nodes;
        for (Eigen::Index node_id = 0; node_id < NumberOfNodesAtCompileTime; ++node_id) {
            nodes.row(node_id) = get_node(node_id);
        }
        return nodes;
    }

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
