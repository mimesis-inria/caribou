#pragma once

#include <Caribou/config.h>
#include <Caribou/Geometry/Base/BaseTetrahedron.h>
#include <Caribou/Geometry/Triangle.h>
#include <Eigen/Core>
#include <array>

namespace caribou::geometry {

struct Tetrahedron;

template<>
struct traits<Tetrahedron> {
    static constexpr UNSIGNED_INTEGER_TYPE CanonicalDimension = 3;
    static constexpr UNSIGNED_INTEGER_TYPE Dimension = 3;
    static constexpr INTEGER_TYPE NumberOfNodesAtCompileTime = 4;
    static constexpr INTEGER_TYPE NumberOfGaussNodesAtCompileTime = 1;

    using BoundaryElementType = Triangle<3>;
    static constexpr INTEGER_TYPE NumberOfBoundaryElementsAtCompileTime = 4;
};

/**
 * Linear Tetrahedron
 *
 * \verbatim
 *                    v
 *                  .
 *                ,/
 *               /
 *            2
 *          ,/|`\
 *        ,/  |  `\
 *      ,/    '.   `\
 *    ,/       |     `\
 *  ,/         |       `\
 * 0-----------'.--------1 --> u
 *  `\.         |      ,/
 *     `\.      |    ,/
 *        `\.   '. ,/
 *           `\. |/
 *              `3
 *                 `\.
 *                    ` w
 * \endverbatim
 *
 */
struct Tetrahedron : public BaseTetrahedron<Tetrahedron> {
    // Types
    using Base = BaseTetrahedron<Tetrahedron>;
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
    Tetrahedron() : Base() {
        this->p_nodes.row(0) = WorldCoordinates(0, 0, 0);
        this->p_nodes.row(1) = WorldCoordinates(1, 0, 0);
        this->p_nodes.row(2) = WorldCoordinates(0, 1, 0);
        this->p_nodes.row(3) = WorldCoordinates(0, 0, 1);
    };

    // Public methods
    /**
     * Get the list of node indices of the edges.
     * \sa Element::boundary_elements_node_indices
     */
    inline auto edges() const {
        static const std::array<std::array<UNSIGNED_INTEGER_TYPE, 2>, 6> indices = {{
             {0, 1}, // Edge 0
             {1, 2}, // Edge 1
             {2, 0}, // Edge 2
             {0, 3}, // Edge 3
             {3, 1}, // Edge 4
             {3, 2}  // Edge 5
         }};
        return indices;
    }


private:
    // Implementations
    friend struct Element<Tetrahedron>;
    friend struct BaseTetrahedron<Tetrahedron>;
    inline auto get_L(const LocalCoordinates & xi) const -> Vector<NumberOfNodesAtCompileTime> {
        const auto  & u = xi[0];
        const auto  & v = xi[1];
        const auto  & w = xi[2];
        return {
            static_cast<FLOATING_POINT_TYPE>(1 - u - v - w),
            static_cast<FLOATING_POINT_TYPE>(u),
            static_cast<FLOATING_POINT_TYPE>(v),
            static_cast<FLOATING_POINT_TYPE>(w)
        };
    };

    inline auto get_dL(const LocalCoordinates & /*xi*/) const {
        Matrix<NumberOfNodesAtCompileTime, CanonicalDimension> m;
        //  dL/du  dL/dv  dL/dw
        m << -1,    -1,    -1,   // Node 0
              1,     0,     0,   // Node 1
              0,     1,     0,   // Node 1
              0,     0,     1;   // Node 3
        return m;
    };

    inline auto get_gauss_nodes() const -> const auto & {
        using Weight = FLOATING_POINT_TYPE;
        static const std::vector<GaussNode> gauss_nodes {
            GaussNode {LocalCoordinates(1/4., 1/4., 1/4.), Weight(1/6.)} // Node 0
        };
        return gauss_nodes;
    }

    inline auto get_boundary_elements_nodes() const -> const auto & {
        static const std::array<std::array<UNSIGNED_INTEGER_TYPE, 3>, 4> indices = {{
            {0, 2, 1}, // Face 0
            {0, 1, 3}, // Face 1
            {0, 3, 2}, // Face 2
            {3, 1, 2}  // Face 3
        }};
        return indices;
    }
};

}
