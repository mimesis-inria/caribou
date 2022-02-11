#pragma once

#include <Caribou/config.h>
#include <Caribou/Geometry/Base/BaseRectangularHexahedron.h>
#include <Caribou/Geometry/Hexahedron.h>
#include <Caribou/Geometry/RectangularQuad.h>
#include <Eigen/Core>

namespace caribou::geometry {

struct RectangularHexahedron20;

template<>
struct traits<RectangularHexahedron20> {
    static constexpr UNSIGNED_INTEGER_TYPE CanonicalDimension = 3;
    static constexpr UNSIGNED_INTEGER_TYPE Dimension = 3;
    static constexpr INTEGER_TYPE NumberOfNodesAtCompileTime = 20;
    static constexpr INTEGER_TYPE NumberOfGaussNodesAtCompileTime = 4;

    using BoundaryElementType = RectangularQuad8<3>;
    static constexpr INTEGER_TYPE NumberOfBoundaryElementsAtCompileTime = 6;
};

/**
 * Quadratic Rectangular Hexahedron (20 nodes)
 *
 * \verbatim
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
 * \endverbatim
 *
 */
struct RectangularHexahedron20 : public BaseRectangularHexahedron<RectangularHexahedron20> {
    // Types
    using Base = BaseRectangularHexahedron<RectangularHexahedron20>;
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

    static constexpr auto canonical_nodes = Hexahedron20::canonical_nodes;

    // Constructors
    using Base::Base;

    // Constructor from a regular Hexahedron
    explicit RectangularHexahedron20(const Hexahedron & hexa) {
        this->p_center = hexa.center();
        this->p_H = Size((hexa.node(1)-hexa.node(0)).norm(), (hexa.node(3)-hexa.node(0)).norm(), (hexa.node(4)-hexa.node(0)).norm());
        this->p_R = hexa.frame({0, 0, 0});
    }

    explicit RectangularHexahedron20(const Hexahedron20 & hexa) {
        this->p_center = hexa.center();
        this->p_H = Size((hexa.node(1)-hexa.node(0)).norm(), (hexa.node(3)-hexa.node(0)).norm(), (hexa.node(4)-hexa.node(0)).norm());
        this->p_R = hexa.frame({0, 0, 0});
    }

    /** Constructor from an Eigen matrix containing the positions of an quadratic hexa nodes */
    template<typename EigenType, REQUIRES(EigenType::RowsAtCompileTime == NumberOfNodesAtCompileTime)>
    explicit RectangularHexahedron20(Eigen::EigenBase<EigenType> & nodes)
    : RectangularHexahedron20(Hexahedron20(nodes)) {}

    /** Constructor from an Eigen matrix containing the positions of an quadratic hexa nodes */
    template<typename EigenType, REQUIRES(EigenType::RowsAtCompileTime == NumberOfNodesAtCompileTime)>
    explicit RectangularHexahedron20(const Eigen::EigenBase<EigenType> & nodes)
    : RectangularHexahedron20(Hexahedron20(nodes)) {}

    /** Constructor from an Eigen matrix containing the positions of a linear hexa nodes */
    template<typename EigenType, REQUIRES(EigenType::RowsAtCompileTime == 8)>
    explicit RectangularHexahedron20(Eigen::EigenBase<EigenType> & nodes)
    : RectangularHexahedron20(Hexahedron(nodes)) {}

    /** Constructor from an Eigen matrix containing the positions of a linear hexa nodes */
    template<typename EigenType, REQUIRES(EigenType::RowsAtCompileTime == 8)>
    explicit RectangularHexahedron20(const Eigen::EigenBase<EigenType> & nodes)
    : RectangularHexahedron20(Hexahedron(nodes)) {}

    // Public methods
    /**
     * Get the list of node indices of the edges.
     * \sa Element::boundary_elements_node_indices
     */
    inline auto edges() const {
        static const std::array<std::array<UNSIGNED_INTEGER_TYPE, 3>, 12> indices = {{
            {0, 1, 8},  // Edge 0
            {1, 2, 9},  // Edge 1
            {2, 3, 10}, // Edge 2
            {3, 0, 11}, // Edge 3
            {4, 5, 12}, // Edge 4
            {5, 6, 13}, // Edge 5
            {6, 7, 14}, // Edge 6
            {7, 4, 15}, // Edge 7
            {0, 4, 16}, // Edge 8
            {1, 5, 17}, // Edge 9
            {2, 6, 18}, // Edge 10
            {3, 7, 19}  // Edge 11
        }};
        return indices;
    }

private:
    // Implementations
    friend struct Element<RectangularHexahedron20>;
    friend struct BaseRectangularHexahedron<RectangularHexahedron20>;
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
