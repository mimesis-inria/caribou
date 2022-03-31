#pragma once

#include <Caribou/config.h>
#include <Caribou/Geometry/Base/BaseHexahedron.h>
#include <Caribou/Geometry/Quad8.h>
#include <Caribou/Geometry/Hexahedron_FEniCS.h>
#include <Eigen/Core>

namespace caribou::geometry {

struct Hexahedron_FEniCS20;

template<>
struct traits<Hexahedron_FEniCS20> {
    static constexpr UNSIGNED_INTEGER_TYPE CanonicalDimension = 3;
    static constexpr UNSIGNED_INTEGER_TYPE Dimension = 3;
    static constexpr INTEGER_TYPE NumberOfNodesAtCompileTime = 20;
    static constexpr INTEGER_TYPE NumberOfGaussNodesAtCompileTime = 8;

    using BoundaryElementType = Quad8<3>;
    static constexpr INTEGER_TYPE NumberOfBoundaryElementsAtCompileTime = 6;
};

/**
 * Quadratic Hexahedron_FEniCS20 (20 nodes)
 *
 * \verbatim
 *         v
 *         ^
 *         |
 *  6-----19---7
 *  |\     |   |\
 *  | 17   |   | 18
 * 14  \   |   15 \
 *  |   4----16+---5
 *  |   |  +-- |-- | ---> u
 *  2---+-13\--3   |
 *   \  10   \  \  12
 *    9 |     \  11|
 *     \|      w  \|
 *      0----8-----1   u
 * \endverbatim
 *
 */
struct Hexahedron_FEniCS20  : public BaseHexahedron<Hexahedron_FEniCS20> {
    // Types
    using Base = BaseHexahedron<Hexahedron_FEniCS20>;
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
        { 0, -1, +1}, // Node 8
        {-1, -1,  0}, // Node 9
        {-1,  0, +1}, // Node 10
        {+1, -1,  0}, // Node 11
        {+1,  0, +1}, // Node 12
        { 0, -1, -1}, // Node 13
        {-1,  0, -1}, // Node 14
        {+1,  0, -1}, // Node 15
        { 0, +1, +1}, // Node 16
        {-1, +1,  0}, // Node 17
        {+1, +1,  0}, // Node 18
        { 0, +1, -1}, // Node 19
    };

    // Constructors
    using Base::Base;
    Hexahedron_FEniCS20() : Base() {
        this->p_nodes.row(0)  = WorldCoordinates(canonical_nodes[0][0],  canonical_nodes[0][1],  canonical_nodes[0][2]);  // Node 0
        this->p_nodes.row(1)  = WorldCoordinates(canonical_nodes[1][0],  canonical_nodes[1][1],  canonical_nodes[1][2]);  // Node 1
        this->p_nodes.row(2)  = WorldCoordinates(canonical_nodes[2][0],  canonical_nodes[2][1],  canonical_nodes[2][2]);  // Node 2
        this->p_nodes.row(3)  = WorldCoordinates(canonical_nodes[3][0],  canonical_nodes[3][1],  canonical_nodes[3][2]);  // Node 3
        this->p_nodes.row(4)  = WorldCoordinates(canonical_nodes[4][0],  canonical_nodes[4][1],  canonical_nodes[4][2]);  // Node 4
        this->p_nodes.row(5)  = WorldCoordinates(canonical_nodes[5][0],  canonical_nodes[5][1],  canonical_nodes[5][2]);  // Node 5
        this->p_nodes.row(6)  = WorldCoordinates(canonical_nodes[6][0],  canonical_nodes[6][1],  canonical_nodes[6][2]);  // Node 6
        this->p_nodes.row(7)  = WorldCoordinates(canonical_nodes[7][0],  canonical_nodes[7][1],  canonical_nodes[7][2]);  // Node 7
        this->p_nodes.row(8)  = WorldCoordinates(canonical_nodes[8][0],  canonical_nodes[8][1],  canonical_nodes[8][2]);  // Node 8
        this->p_nodes.row(9)  = WorldCoordinates(canonical_nodes[9][0],  canonical_nodes[9][1],  canonical_nodes[9][2]);  // Node 9
        this->p_nodes.row(10) = WorldCoordinates(canonical_nodes[10][0], canonical_nodes[10][1], canonical_nodes[10][2]); // Node 10
        this->p_nodes.row(11) = WorldCoordinates(canonical_nodes[11][0], canonical_nodes[11][1], canonical_nodes[11][2]); // Node 11
        this->p_nodes.row(12) = WorldCoordinates(canonical_nodes[12][0], canonical_nodes[12][1], canonical_nodes[12][2]); // Node 12
        this->p_nodes.row(13) = WorldCoordinates(canonical_nodes[13][0], canonical_nodes[13][1], canonical_nodes[13][2]); // Node 13
        this->p_nodes.row(14) = WorldCoordinates(canonical_nodes[14][0], canonical_nodes[14][1], canonical_nodes[14][2]); // Node 14
        this->p_nodes.row(15) = WorldCoordinates(canonical_nodes[15][0], canonical_nodes[15][1], canonical_nodes[15][2]); // Node 15
        this->p_nodes.row(16) = WorldCoordinates(canonical_nodes[16][0], canonical_nodes[16][1], canonical_nodes[16][2]); // Node 16
        this->p_nodes.row(17) = WorldCoordinates(canonical_nodes[17][0], canonical_nodes[17][1], canonical_nodes[17][2]); // Node 17
        this->p_nodes.row(18) = WorldCoordinates(canonical_nodes[18][0], canonical_nodes[18][1], canonical_nodes[18][2]); // Node 18
        this->p_nodes.row(19) = WorldCoordinates(canonical_nodes[19][0], canonical_nodes[19][1], canonical_nodes[19][2]); // Node 19

    };

    // Construct a quadratic Hexahedron_FEniCS from a linear one
    explicit Hexahedron_FEniCS20(const Hexahedron_FEniCS & linear_Hexahedron_FEniCS) : Base() {
        this->p_nodes.row(0) = linear_Hexahedron_FEniCS.node(0); // Node 0
        this->p_nodes.row(1) = linear_Hexahedron_FEniCS.node(1); // Node 1
        this->p_nodes.row(2) = linear_Hexahedron_FEniCS.node(2); // Node 2
        this->p_nodes.row(3) = linear_Hexahedron_FEniCS.node(3); // Node 3
        this->p_nodes.row(4) = linear_Hexahedron_FEniCS.node(4); // Node 4
        this->p_nodes.row(5) = linear_Hexahedron_FEniCS.node(5); // Node 5
        this->p_nodes.row(6) = linear_Hexahedron_FEniCS.node(6); // Node 6
        this->p_nodes.row(7) = linear_Hexahedron_FEniCS.node(7); // Node 7
        this->p_nodes.row(8)  = linear_Hexahedron_FEniCS.world_coordinates(LocalCoordinates(canonical_nodes[8][0],  canonical_nodes[8][1],  canonical_nodes[8][2]));  // Node 8
        this->p_nodes.row(9)  = linear_Hexahedron_FEniCS.world_coordinates(LocalCoordinates(canonical_nodes[9][0],  canonical_nodes[9][1],  canonical_nodes[9][2]));  // Node 9
        this->p_nodes.row(10) = linear_Hexahedron_FEniCS.world_coordinates(LocalCoordinates(canonical_nodes[10][0], canonical_nodes[10][1], canonical_nodes[10][2])); // Node 10
        this->p_nodes.row(11) = linear_Hexahedron_FEniCS.world_coordinates(LocalCoordinates(canonical_nodes[11][0], canonical_nodes[11][1], canonical_nodes[11][2])); // Node 11
        this->p_nodes.row(12) = linear_Hexahedron_FEniCS.world_coordinates(LocalCoordinates(canonical_nodes[12][0], canonical_nodes[12][1], canonical_nodes[12][2])); // Node 12
        this->p_nodes.row(13) = linear_Hexahedron_FEniCS.world_coordinates(LocalCoordinates(canonical_nodes[13][0], canonical_nodes[13][1], canonical_nodes[13][2])); // Node 13
        this->p_nodes.row(14) = linear_Hexahedron_FEniCS.world_coordinates(LocalCoordinates(canonical_nodes[14][0], canonical_nodes[14][1], canonical_nodes[14][2])); // Node 14
        this->p_nodes.row(15) = linear_Hexahedron_FEniCS.world_coordinates(LocalCoordinates(canonical_nodes[15][0], canonical_nodes[15][1], canonical_nodes[15][2])); // Node 15
        this->p_nodes.row(16) = linear_Hexahedron_FEniCS.world_coordinates(LocalCoordinates(canonical_nodes[16][0], canonical_nodes[16][1], canonical_nodes[16][2])); // Node 16
        this->p_nodes.row(17) = linear_Hexahedron_FEniCS.world_coordinates(LocalCoordinates(canonical_nodes[17][0], canonical_nodes[17][1], canonical_nodes[17][2])); // Node 17
        this->p_nodes.row(18) = linear_Hexahedron_FEniCS.world_coordinates(LocalCoordinates(canonical_nodes[18][0], canonical_nodes[18][1], canonical_nodes[18][2])); // Node 18
        this->p_nodes.row(19) = linear_Hexahedron_FEniCS.world_coordinates(LocalCoordinates(canonical_nodes[19][0], canonical_nodes[19][1], canonical_nodes[19][2])); // Node 19
    }

    /** Construct a quadratic Hexahedron_FEniCS20 from eight nodes */
    Hexahedron_FEniCS20(WorldCoordinates & p0, WorldCoordinates & p1, WorldCoordinates & p2, WorldCoordinates & p3, WorldCoordinates & p4, WorldCoordinates & p5, WorldCoordinates & p6, WorldCoordinates & p7)
        : Hexahedron_FEniCS20(Hexahedron_FEniCS(p0, p1, p2, p3, p4, p5, p6, p7)) {}

    /** Construct a quadratic Hexahedron_FEniCS20 from eight nodes */
    Hexahedron_FEniCS20(const WorldCoordinates & p0, const WorldCoordinates & p1, const WorldCoordinates & p2, const WorldCoordinates & p3, const WorldCoordinates & p4, const WorldCoordinates & p5, const WorldCoordinates & p6, const WorldCoordinates & p7)
        : Hexahedron_FEniCS20(Hexahedron_FEniCS(p0, p1, p2, p3, p4, p5, p6, p7)) {}

    // Public methods
    /**
     * Get the list of node indices of the edges.
     * \sa Element::boundary_elements_node_indices
     */
    inline auto edges() const {
        static const std::array<std::array<UNSIGNED_INTEGER_TYPE, 3>, 12> indices = {{
             {2, 13, 3},  // Edge 0
             {3, 7, 15},  // Edge 1
             {7, 6, 19}, // Edge 2
             {6, 2, 14}, // Edge 3
             {0, 1,  8}, // Edge 4
             {1, 5, 12}, // Edge 5
             {5, 4, 16}, // Edge 6
             {4, 0, 10}, // Edge 7
             {2, 0,  9}, // Edge 8
             {3, 1, 11}, // Edge 9
             {7, 5, 18}, // Edge 10
             {6, 4, 17}  // Edge 11
         }};
        return indices;
    }

private:
    // Implementations
    friend struct Element<Hexahedron_FEniCS20>;
    friend struct BaseHexahedron<Hexahedron_FEniCS20>;
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
            {2, 6, 7, 3, 14, 19, 15, 13}, // Face 0
            {2, 0, 4, 6,  9, 10, 17, 14}, // Face 1
            {3, 7, 5, 1, 15, 18, 12, 11}, // Face 2
            {2, 3, 1, 0, 13, 11,  8,  9}, // Face 3
            {7, 6, 4, 5, 19, 17, 16, 18}, // Face 4
            {0, 1, 5, 4,  8, 12, 16, 10}  // Face 5
        }};
        return indices;
    }
};

}
