#pragma once

#include <Caribou/config.h>
#include <Caribou/Geometry/Base/BaseTetrahedron.h>
#include <Caribou/Geometry/Triangle6.h>
#include <Caribou/Geometry/Tetrahedron.h>
#include <Eigen/Core>
#include <array>

namespace caribou::geometry {

struct Tetrahedron10;

template<>
struct traits<Tetrahedron10> {
    static constexpr UNSIGNED_INTEGER_TYPE CanonicalDimension = 3;
    static constexpr UNSIGNED_INTEGER_TYPE Dimension = 3;
    static constexpr INTEGER_TYPE NumberOfNodesAtCompileTime = 10;
    static constexpr INTEGER_TYPE NumberOfGaussNodesAtCompileTime = 4;

    using BoundaryElementType = Triangle6<3>;
    static constexpr INTEGER_TYPE NumberOfBoundaryElementsAtCompileTime = 4;
};

/**
 * Quadratic Tetrahedron
 *
 * \verbatim
 *                    v
 *                  .
 *                ,/
 *               /
 *            2
 *          ,/|`\
 *        ,/  |  `\
 *      ,6    '.   `5
 *    ,/       9     `\
 *  ,/         |       `\
 * 0--------4--'.--------1 --> u
 *  `\.         |      ,/
 *     `\.      |    ,8
 *        `7.   '. ,/
 *           `\. |/
 *              `3
 *                 `\.
 *                    ` w
 * \endverbatim
 *
 */
struct Tetrahedron10: public BaseTetrahedron<Tetrahedron10> {
    // Types
    using Base = BaseTetrahedron<Tetrahedron10>;
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
    Tetrahedron10() : Base() {
        this->p_nodes.row(0) = WorldCoordinates(0.0, 0.0, 0.0); // Node 0
        this->p_nodes.row(1) = WorldCoordinates(1.0, 0.0, 0.0); // Node 1
        this->p_nodes.row(2) = WorldCoordinates(0.0, 1.0, 0.0); // Node 2
        this->p_nodes.row(3) = WorldCoordinates(0.0, 0.0, 1.0); // Node 3
        this->p_nodes.row(4) = WorldCoordinates(0.5, 0.0, 0.0); // Node 4
        this->p_nodes.row(5) = WorldCoordinates(0.5, 0.5, 0.0); // Node 5
        this->p_nodes.row(6) = WorldCoordinates(0.0, 0.5, 0.0); // Node 6
        this->p_nodes.row(7) = WorldCoordinates(0.0, 0.0, 0.5); // Node 7
        this->p_nodes.row(8) = WorldCoordinates(0.5, 0.0, 0.5); // Node 8
        this->p_nodes.row(9) = WorldCoordinates(0.0, 0.5, 0.5); // Node 9
    };

    // Construct a quadratic Tetrahedron from a linear one
    explicit Tetrahedron10(const Tetrahedron & linear_Tetrahedron) : Base() {
        this->p_nodes.row(0) = linear_Tetrahedron.node(0); // Node 0
        this->p_nodes.row(1) = linear_Tetrahedron.node(1); // Node 1
        this->p_nodes.row(2) = linear_Tetrahedron.node(2); // Node 2
        this->p_nodes.row(3) = linear_Tetrahedron.node(3); // Node 3
        this->p_nodes.row(4) = linear_Tetrahedron.world_coordinates(LocalCoordinates(0.5, 0.0, 0.0)); // Node 4
        this->p_nodes.row(5) = linear_Tetrahedron.world_coordinates(LocalCoordinates(0.5, 0.5, 0.0)); // Node 5
        this->p_nodes.row(6) = linear_Tetrahedron.world_coordinates(LocalCoordinates(0.0, 0.5, 0.0)); // Node 6
        this->p_nodes.row(7) = linear_Tetrahedron.world_coordinates(LocalCoordinates(0.0, 0.0, 0.5)); // Node 7
        this->p_nodes.row(8) = linear_Tetrahedron.world_coordinates(LocalCoordinates(0.5, 0.0, 0.5)); // Node 8
        this->p_nodes.row(9) = linear_Tetrahedron.world_coordinates(LocalCoordinates(0.0, 0.5, 0.5)); // Node 9
    }

    /** Construct a quadratic Tetrahedron from four nodes */
    Tetrahedron10(WorldCoordinates & p0, WorldCoordinates & p1, WorldCoordinates & p2, WorldCoordinates & p3)
    : Tetrahedron10(Tetrahedron(p0, p1, p2, p3)) {}

    /** Construct a quadratic Tetrahedron from four nodes */
    Tetrahedron10(const WorldCoordinates & p0, const WorldCoordinates & p1, const WorldCoordinates & p2, const WorldCoordinates & p3)
    : Tetrahedron10(Tetrahedron(p0, p1, p2, p3)) {}

    // Public methods
    /**
     * Get the list of node indices of the edges.
     * \sa Element::boundary_elements_node_indices
     */
    inline auto edges() const {
        static const std::array<std::array<UNSIGNED_INTEGER_TYPE, 3>, 6> indices = {{
            {0, 1, 4}, // Edge 0
            {1, 2, 5}, // Edge 1
            {2, 0, 6}, // Edge 2
            {0, 3, 7}, // Edge 3
            {3, 1, 8}, // Edge 4
            {3, 2, 9}  // Edge 5
        }};
        return indices;
    }

private:
    // Implementations
    friend struct Element<Tetrahedron10>;
    friend struct BaseTetrahedron<Tetrahedron10>;
    inline auto get_L(const LocalCoordinates & xi) const {
        const auto  & u = xi[0];
        const auto  & v = xi[1];
        const auto  & w = xi[2];
        const auto l = 1 - u - v - w;

        Vector<NumberOfNodesAtCompileTime> m;
        m << l * (2*l - 1),   // Node 0
             u * (2*u - 1),   // Node 1
             v * (2*v - 1),   // Node 2
             w * (2*w - 1),   // Node 3
             4 * l * u,       // Node 4
             4 * u * v,       // Node 5
             4 * v * l,       // Node 6
             4 * l * w,       // Node 7
             4 * u * w,       // Node 8
             4 * v * w;       // Node 9

        return m;
    };

    inline auto get_dL(const LocalCoordinates & xi) const {
        const auto & u = xi[0];
        const auto & v = xi[1];
        const auto  & w = xi[2];
        const auto l = 1 - u - v - w;

        Matrix<NumberOfNodesAtCompileTime, CanonicalDimension> m;
        //       dL/du              dL/dv              dL/dw
        m <<   1 - (4 * l),       1 - (4 * l),       1 - (4 * l),   // Node 0
              (4 * u) - 1 ,          0       ,          0       ,   // Node 1
                    0     ,      (4 * v) - 1 ,          0       ,   // Node 2
                    0     ,          0       ,      (4 * w) - 1 ,   // Node 3
               4 * (l - u),      -4 * u      ,      -4 * u      ,   // Node 4
                4 * v     ,       4 * u      ,          0       ,   // Node 5
               -4 * v     ,       4 * (l - v),      -4 * v      ,   // Node 6
               -4 * w     ,      -4 * w      ,       4 * (l - w),   // Node 7
                4 * w     ,          0       ,       4 * u      ,   // Node 8
                    0     ,       4 * w      ,       4 * v      ;   // Node 9
        return m;
    };

    inline auto get_gauss_nodes() const -> const auto & {
        using Weight = FLOATING_POINT_TYPE;
        static const std::vector<GaussNode> gauss_nodes {
            GaussNode {LocalCoordinates(0.1381966011250105, 0.1381966011250105, 0.1381966011250105), Weight(1/24.)}, // Node 0
            GaussNode {LocalCoordinates(0.1381966011250105, 0.1381966011250105, 0.5854101966249685), Weight(1/24.)}, // Node 1
            GaussNode {LocalCoordinates(0.1381966011250105, 0.5854101966249685, 0.1381966011250105), Weight(1/24.)}, // Node 2
            GaussNode {LocalCoordinates(0.5854101966249685, 0.1381966011250105, 0.1381966011250105), Weight(1/24.)}  // Node 3
        };
        return gauss_nodes;
    }

    inline auto get_boundary_elements_nodes() const -> const auto & {
        static const std::array<std::array<UNSIGNED_INTEGER_TYPE, 6>, 4> indices = {{
            {0, 2, 1, 6, 5, 4}, // Face 0
            {0, 1, 3, 4, 8, 7}, // Face 1
            {0, 3, 2, 7, 9, 6}, // Face 2
            {3, 1, 2, 8, 5, 9}  // Face 3
        }};
        return indices;
    }
};

}
