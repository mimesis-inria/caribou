#ifndef CARIBOU_TOPOLOGY_GRID_INTERNAL_GRID_H
#define CARIBOU_TOPOLOGY_GRID_INTERNAL_GRID_H

#include <cstddef>
#include <list>
#include <array>
#include <vector>
#include <bitset>
#include <Caribou/config.h>
#include <Caribou/macros.h>
#include <Eigen/Core>

namespace caribou::topology::internal {

/**
 * Simple representation of a Grid in space.
 *
 * ** Do not use this class directly. Use instead caribou::topology::Grid. **
 *
 * The functions declared in this class can be used with any type of grids (static grid, container grid, etc.).
 *
 * Do to so, it uses the Curiously recurring template pattern (CRTP) :
 *    https://en.wikipedia.org/wiki/Curiously_recurring_template_pattern
 *
 * A Grid is a set of multiple cell entities (same-length lines in 1D, rectangles
 * in 2D and rectangular hexahedrons in 3D) aligned in the x, y and z axis.
 *
 * @tparam Dim Dimension of the grid (1D, 2D or 3D).
 * @tparam GridType_ Type of the derived grid class that will implement the final functions.
 */
template <size_t Dim, class GridType_>
struct BaseGrid
{
    static constexpr size_t Dimension = Dim;

    using GridType = GridType_;

    using Int = INTEGER_TYPE;
    using UInt = UNSIGNED_INTEGER_TYPE;
    using Float = FLOATING_POINT_TYPE;

    using VecFloat = Eigen::Matrix<Float, Dimension, 1>;
    using VecInt = Eigen::Matrix<Int, Dimension, 1>;
    using VecUInt = Eigen::Matrix<UInt, Dimension, 1>;

    using NodeIndex = Int;
    using EdgeIndex = Int;
    using CellIndex = Int;
    using Dimensions = VecFloat;
    using Subdivisions = VecUInt;
    using LocalCoordinates = VecFloat;
    using WorldCoordinates = VecFloat;
    using GridCoordinates = VecInt;
    using CellSet = std::list<CellIndex>;

    static_assert(Dimension == 1 || Dimension == 2 || Dimension == 3, "Grids are only available in 1, 2 or 3 dimensions");


    /** Default constructor **/
    constexpr
    BaseGrid() = delete;

    /**
     * Main constructor of a grid
     * @param anchor_position Position of the anchor point. The anchor point is the corner node indexed with 0
     * when seeing the grid as a big regular hexahedral element.
     * @param subdivisions Number of sub-cells in each directions of the grid.
     * @param dimension Global size of the grid.
     */
    constexpr
    BaseGrid(const WorldCoordinates & anchor_position, const Subdivisions & n, const Dimensions & size)
            : m_anchor_position(anchor_position), m_number_of_subdivisions(n), m_size(size)
    {}

    /** Get the number of cells (same-length lines in 1D, rectangles
     * in 2D and rectangular hexahedrons in 3D) in this grid. **/
    [[nodiscard]] inline auto
    number_of_cells() const noexcept -> UInt
    {
        const auto & n = m_number_of_subdivisions;
        return n.prod();
    }

    /** Get the number of distinct nodes in this grid. **/
    [[nodiscard]] inline auto
    number_of_nodes() const noexcept -> UInt
    {
        const auto & n = m_number_of_subdivisions;
        const auto & one = VecUInt::Ones();

        return (n+one).prod();
    }

    /** Get position of the anchor point of this grid. The anchor point is the corner node indexed with 0
     * when seeing the grid as a big regular hexahedral element. */
    inline auto
    anchor_position () const noexcept -> const WorldCoordinates &
    {
        return m_anchor_position;
    }

    /** Get the number of sub-cells in this grid : nx (in 2D), {nx, ny} (in 2D) or {nx, ny, nz} (in 3D) */
    inline auto
    N () const noexcept -> const Subdivisions &
    {
        return m_number_of_subdivisions;
    }

    /** Get the dimension of a cell : hx (in 1D), {hx, hy} (in 2D) or {hx, hy, hz} (in 3D) */
    inline auto
    H () const noexcept -> Dimensions
    {
        return (m_size.array() / m_number_of_subdivisions.array().template cast<FLOATING_POINT_TYPE>()).matrix();
    }

    /** The global dimension of this grid : hx (in 1D), {hx, hy} (in 2D) or {hx, hy, hz} (in 3D) */
    inline auto
    size() const noexcept -> Dimensions
    {
        return m_size;
    }

    /** Get the index of the cell that contains the given world coordinates.
     *
     * @param test_for_boundary If the location of the point lies directly on one of the grid's upper boundary
     * (p_x == grid_max_x or p_y == grid_max_y or p_z == grid_max_z), the return grid coordinate will be
     * outside of the grid boundaries. Enabling this option will verify that and return the closest cell inside.
     * */
    inline auto
    cell_index_containing(const WorldCoordinates & coordinates, bool test_for_boundary = false) const noexcept -> CellIndex
    {
        return Self().cell_index_at(cell_coordinates_at(coordinates, test_for_boundary));
    }

    /**
     * Get all the cells around the given world coordinates.
     * If the position is inside a cell (1D, 2D or 3D), only one cell is returned.
     * If it is on a node (1D, 2D or 3D), 2 cells are returned in 1D, 4 cells in 2D, and 8 cells in 3D.
     * If it is on an edge (2D and 3D), two cells are returned in 2D, four in 3D.
     * If it is on a face (3D), two cells are returned in 3D.
     */
    inline auto
    cells_around(const WorldCoordinates & coordinates) const noexcept -> std::vector<CellIndex>
    {

        struct CwiseRound {
             auto operator()(const Float& x) const -> Float { return std::round(x); }
        };
        struct CwiseFloor {
            auto operator()(const Float& x) const -> Float { return std::floor(x); }
        };


        std::vector<CellIndex> cells;
        cells.reserve((unsigned) 1<<Dimension);

        const CellIndex ncells = number_of_cells();
        const VecFloat h = H();
        const VecFloat absolute = ((coordinates - m_anchor_position).array() / h.array()).matrix();
        const VecFloat rounded = absolute.unaryExpr(CwiseRound());
        const VecFloat distance = absolute - rounded;

        auto close_to_axis = std::bitset<Dimension>();

        // Let's find the axis for which our coordinate is very close to
        for (UNSIGNED_INTEGER_TYPE axis = 0; axis < Dimension; ++axis) {
            if (distance[axis]*distance[axis] < EPSILON*EPSILON )
                close_to_axis[axis] = true;
        }

        if (close_to_axis.none()) {
            // We are not near any axis, which means we are well inside a cell's boundaries
            const auto cell_index = Self().cell_index_at(absolute.unaryExpr(CwiseFloor()). template cast<Int>());
            if (IN_CLOSED_INTERVAL(0, cell_index, ncells-1)) {
                cells.emplace_back(cell_index);
            }
        } else {
            const VecFloat d = VecFloat::Constant(0.5);
            const auto & n = N();

            std::array<std::vector<UNSIGNED_INTEGER_TYPE>, Dimension> axis_indices;
            for (auto & indices  : axis_indices) {
                indices.reserve(Dimension);
            }

            for (UNSIGNED_INTEGER_TYPE axis = 0; axis < Dimension; ++axis) {
                if (close_to_axis[axis]) {
                    auto index = floor(rounded[axis]-d[axis]);
                    if (index >= 0) {
                        axis_indices[axis].emplace_back(index);
                    }
                    index = floor(rounded[axis]+d[axis]);
                    if (index <= n[axis]-1) {
                        axis_indices[axis].emplace_back(index);
                    }
                } else {
                    auto index = floor(absolute[axis]);
                    if (IN_CLOSED_INTERVAL(0, index, n[axis]-1)) {
                        axis_indices[axis].emplace_back(index);
                    }
                }
            }

            for (Eigen::Index i : axis_indices[0]) {
                if constexpr (Dimension == 1) {
                    cells.emplace_back(Self().cell_index_at(GridCoordinates(i)));
                } else {
                    for (Eigen::Index j : axis_indices[1]) {
                        if constexpr (Dimension == 2) {
                            cells.emplace_back(Self().cell_index_at({i, j}));
                        } else {
                            for (Eigen::Index k : axis_indices[2]) {
                                cells.emplace_back(Self().cell_index_at({i, j, k}));
                            }
                        }
                    }
                }
            }
        }

        return cells;
    }

    /** Get the grid location of the cell that contains the given world coordinates
     *
     * @param test_for_boundary If the location of the point lies directly on one of the grid's upper boundary
     * (p_x == grid_max_x or p_y == grid_max_y or p_z == grid_max_z), the return grid coordinate will be
     * outside of the grid boundaries. Enabling this option will verify that and return the closest cell inside.
     * */
    inline auto
    cell_coordinates_at(const WorldCoordinates & coordinates, bool test_for_boundary) const noexcept -> GridCoordinates
    {
        // @todo (jnbrunet2000@gmail.com): This test should be using inverse mapping function from a regular
        //  hexahedron geometric element defined in the geometry module.

        if (test_for_boundary) {
            WorldCoordinates p = coordinates;
            const auto h = H();
            const WorldCoordinates epsilon = h/1000000;
            const WorldCoordinates margin = h/10.;
            const WorldCoordinates upper_bound = m_anchor_position + (N().array(). template cast<FLOATING_POINT_TYPE>()*h.array()).matrix();
            for (UNSIGNED_INTEGER_TYPE i = 0; i < Dimension; ++i) {
                if (IN_CLOSED_INTERVAL(-epsilon[i], coordinates[i] - upper_bound[i], epsilon[i])) {
                    p[i] -= margin[i];
                }
            }
            return ((p - m_anchor_position).array() / h.array()).matrix(). template cast<Int>();
        }

        return ((coordinates - m_anchor_position).array() / H().array()).matrix(). template cast<Int>();
    }

    /** Get the grid location of the cell that contains the given world coordinates */
    inline auto
    cell_coordinates_at(const WorldCoordinates & coordinates) const noexcept -> GridCoordinates
    {
        // @todo (jnbrunet2000@gmail.com): This test should be using inverse mapping function from a regular
        //  hexahedron geometric element defined in the geometry module.

        return ((coordinates - m_anchor_position).array() / H().array()).matrix(). template cast<Int>();
    }

    /** Get the node location in world coordinates at grid coordinates */
    inline auto
    node(const GridCoordinates & coordinates) const noexcept -> WorldCoordinates
    {
        const WorldCoordinates relative_position = coordinates. template cast<FLOATING_POINT_TYPE>().cwiseProduct(H());
        return m_anchor_position + relative_position;
    }

    /** Test if the position of a given point is inside this grid */
    inline auto
    contains(const WorldCoordinates & coordinates, const Float epsilon = EPSILON) const noexcept -> bool
    {
        // @todo (jnbrunet2000@gmail.com): This test should be using inverse mapping function from a regular
        //  hexahedron geometric element defined in the geometry module.

        const VecFloat distance_to_anchor = ((coordinates - m_anchor_position).array() / size().array()).matrix();

        if (distance_to_anchor[0] < 0 - epsilon || distance_to_anchor[0] > 1 + epsilon)
            return false;

        if CONSTEXPR_IF (Dimension >= 2) {
            if (distance_to_anchor[1] < 0 - epsilon || distance_to_anchor[1] > 1 + epsilon)
                return false;
        }

        if CONSTEXPR_IF (Dimension == 3) {
            if (distance_to_anchor[2] < 0 - epsilon || distance_to_anchor[2] > 1 + epsilon) {
                return false;
            }
        }

        return true;
    }

    /**
     * Returns the set of cells that enclose (in a bounding-box manner) the given world positions.
     *
     * If one or more positions are found outside the grid, the clipped bounding box will be returned: the bounding box
     * cells will stop at the grid boundaries (no invalid cells will be returned). If all positions are found outside the
     * grid, an empty set is returned.
     *
     * */
    template<typename ...WorldCoordinatesTypes>
    inline auto
    cells_enclosing(const WorldCoordinates & first_position, WorldCoordinatesTypes && ... remaining_positions) const noexcept -> CellSet
    {
        std::array<WorldCoordinates , sizeof...(remaining_positions)+1> positions {{
            first_position, std::forward<WorldCoordinates>(remaining_positions)...
        }};

        // Grid's boundaries (the enclosing cells must not exceed these boundaries)
        GridCoordinates lower_grid_boundary = GridCoordinates::Zero();
        GridCoordinates upper_grid_boundary = N(). template cast<Int>() - GridCoordinates::Ones();

        // Tool function that returns the minimum grid coordinates between two grid coordinates
        auto min_between =
                [this, &lower_grid_boundary, &upper_grid_boundary]
                (const std::vector<CellIndex> & cells_indices, const GridCoordinates & min_coordinates) -> GridCoordinates
                {
                    GridCoordinates lower_bound = min_coordinates;

                    for (const auto & cell_index: cells_indices) {
                        const GridCoordinates coordinates = Self().cell_coordinates_at(cell_index);
                        for (size_t i = 0; i < Dimension; ++i) {
                            lower_bound[i] = std::min(coordinates[i], lower_bound[i]);

                            // Clip to grid boundaries if the coordinates lies outside of the grid
                            lower_bound[i] = std::max(lower_bound[i], lower_grid_boundary[i]);
                            lower_bound[i] = std::min(lower_bound[i], upper_grid_boundary[i]);
                        }
                    }
                    return lower_bound;
                };

        // Tool function that returns the maximum grid coordinates between two grid coordinates
        auto max_between =
                [this, &lower_grid_boundary, &upper_grid_boundary]
                (const std::vector<CellIndex> & cells_indices, const GridCoordinates & max_coordinates) -> GridCoordinates
                {
                    GridCoordinates upper_bound = max_coordinates;

                    for (const auto & cell_index: cells_indices) {
                        const GridCoordinates coordinates = Self().cell_coordinates_at(cell_index);
                        for (size_t i = 0; i < Dimension; ++i) {

                            upper_bound[i] = std::max(coordinates[i], upper_bound[i]);

                            // Clip to grid boundaries if the coordinate lies outside of the grid
                            upper_bound[i] = std::max(upper_bound[i], lower_grid_boundary[i]);
                            upper_bound[i] = std::min(upper_bound[i], upper_grid_boundary[i]);
                        }
                    }
                    return upper_bound;
                };

        // 1. Find the grid coordinates bounding box of the cells that contain each nodes
        GridCoordinates lower_cell = min_between(cells_around(positions[0]), upper_grid_boundary);
        GridCoordinates upper_cell = max_between(cells_around(positions[0]), lower_grid_boundary);

        // This is used to discard the bbox computation if all points are outside the grid
        bool grid_contains_at_least_one_point = contains(WorldCoordinates(positions[0]));

        for (size_t i = 1; i < positions.size(); ++ i) {

            if (not grid_contains_at_least_one_point and contains(WorldCoordinates(positions[i])))
                grid_contains_at_least_one_point = true;

            const auto coordinates = cells_around(positions[i]);

            lower_cell = min_between(coordinates, lower_cell);
            upper_cell = max_between(coordinates, upper_cell);
        }

        if (not grid_contains_at_least_one_point)
            return {};

        // 2. Append all cells within the bounding-box
        CellSet enclosing_cells;
        if CONSTEXPR_IF (Dimension ==1) {
            for (CellIndex i = lower_cell[0]; i <= upper_cell[0]; ++i)
                enclosing_cells.emplace_back(Self().cell_index_at(GridCoordinates(i)));
        } else if CONSTEXPR_IF (Dimension == 2) {
            for (CellIndex j = lower_cell[1]; j <= upper_cell[1]; ++j)
                for (CellIndex i = lower_cell[0]; i <= upper_cell[0]; ++i)
                    enclosing_cells.emplace_back(Self().cell_index_at({i, j}));
        } else { // Dimension == 3
            for (CellIndex k = lower_cell[2]; k <= upper_cell[2]; ++k)
                for (CellIndex j = lower_cell[1]; j <= upper_cell[1]; ++j)
                    for (CellIndex i = lower_cell[0]; i <= upper_cell[0]; ++i)
                        enclosing_cells.emplace_back(Self().cell_index_at({i, j, k}));
        }

        return enclosing_cells;
    }

protected:
    ///< Position of the anchor point of this grid. The anchor point is the corner node indexed with 0
    ///< when seeing the grid as a big regular hexahedral element.
    const WorldCoordinates m_anchor_position;

    ///< Number of sub-cells in the x, y and z directions respectively.
    const Subdivisions m_number_of_subdivisions;

    ///<  Dimension of the grid from the anchor point in the x, y and z directions respectively.
    const Dimensions m_size;

private:
    const inline GridType &
    Self() const
    {
        return static_cast<const GridType &> (*this);
    }

};

} // namespace caribou::topology::internal

#endif //CARIBOU_TOPOLOGY_GRID_INTERNAL_GRID_H
