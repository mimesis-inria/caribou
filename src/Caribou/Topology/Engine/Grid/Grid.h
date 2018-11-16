#ifndef CARIBOU_TOPOLOGY_ENGINE_GRID_GRID_H
#define CARIBOU_TOPOLOGY_ENGINE_GRID_GRID_H

#include <vector>
#include <memory>
#include <cmath>

#include <Caribou/config.h>
#include <Caribou/Geometry/Hexahedron.h>
#include <Caribou/Geometry/Quad.h>


namespace caribou
{

namespace topology
{

namespace engine
{

/**
 * A Grid is a rectangular 2D quad (resp. 3D hexahedron) that contain multiple Cell entities aligned in the x, y (and z in 3D) axis.
 * @tparam TCell The type of cell this grid contains.
 */
template <class TCell>
struct Grid
{
    using CellType = TCell;
    using VecFloat = typename CellType::VecFloat;
    using Index = typename CellType::Index;
    using VecInt = typename CellType::VecInt;
    using Int = typename VecInt::ValueType;
    using Float = typename VecFloat::ValueType;

    static constexpr size_t Dimension = CellType::Dimension;
    static constexpr size_t NumberOfNodes = CellType::NumberOfNodes;

    struct CellRange;

    /** Default constructor is not permitted **/
    Grid() = delete;

    /**
     * Constructor of a grid
     * @param anchor The anchor point position vector (x, y, z). It should be positioned at the center of the grid.
     * @param subdivisions Vector of 3 integer (nx, ny, nz) which specify the number of sub-cells in the x, y and z directions respectively
     * @param dimensions Vector of 3 float (sx, sy, sz) which specify the dimension of the grid from the anchor point in the x, y and z directions respectively
     * @param parent If null, the grid will be initialized as the top level grid. Else, the grid is a subcell of the Grid parent.
     */
    Grid(VecFloat anchor, VecInt subdivisions, VecFloat dimensions);

    /** Get the number of cell subdivisions (nx, ny, nz) of this grid. **/
    inline VecInt
    number_of_subdivision() const
    {return m_number_of_subdivisions;};

    /** Get this grid dimensions (sx, sy, sz). **/
    inline VecFloat
    size() const
    { return m_dimensions;};

    /** Get dimensions (hx, hy, hz) of a top-cell in this grid. **/
    inline VecFloat
    cell_size() const
    {
        const auto & N = number_of_subdivision();
        const auto & S = size();
        const auto   H = S.direct_division(N);

        return H;
    };

    /** Get dimensions (hx, hy, hz) of a given (sub or top) cell in this grid. **/
    inline VecFloat
    cell_size(const CellType & cell) const
    {
        const auto H = cell_size(); // Top-cell size
        const auto l = cell.level(); // Cell level
        return H / std::pow(2, cell.level());
    };

    /**
     * Get the cell located at grid_index (i, j, k).
     * @param grid_coordinates Cell location provided in terms of grid coordinates (i, j, k).
     * @throws std::out_of_range when the grid coordinates are outside of this grid subdivisions (nx, ny, nz).
     */
    inline CellType &
    get(const VecInt & grid_coordinates)
    {
        return m_cells[cell_index(grid_coordinates)];
    }

    /**
     * Get the cell located at grid_index (i, j, k).
     * @param grid_coordinates Cell location provided in terms of grid coordinates (i, j, k).
     * @throws std::out_of_range when the grid coordinates are outside of this grid subdivisions (nx, ny, nz).
     */
    inline const CellType &
    get(const VecInt & grid_coordinates) const
    {
        return get(grid_coordinates);
    }

    /**
     * Get the cell located at a given index.
     * @throws std::out_of_range when the index is outside of this grid range of cells.
     */
    inline CellType &
    get(const Int & index)
    {
        if (index >= m_cells.size())
            throw std::out_of_range(
                    "Trying to access a cell at index " + std::to_string(index) +
                    " where this grid contains only " << std::to_string(m_cells.size()) << " outer cells.");
        return m_cells[index];
    }

    /**
     * Get the cell located at a given index.
     * @throws std::out_of_range when the index is outside of this grid range of cells.
     */
    inline const CellType &
    get(const Int & index) const
    {
        return get(index);
    }

    /**
     * Get the cell index at grid location (i, j, k).
     * @param grid_coordinates Cell location provided in terms of grid coordinates (i, j, k).
     * @throws std::out_of_range when the grid coordinates are outside of this grid subdivisions (nx, ny, nz).
     */
    virtual Index cell_index(const VecInt & grid_coordinates) const;

    /** Get the grid location (i, j, k) at index cell_index */
    virtual VecInt grid_coordinates(const Index & cell_index) const;

    /**
     * Get the indices of the nodes forming the cell located at grid_index (i, j, k). The indices of the nodes will be return
     * following caribou::geometry::LinearHexahedron::node node numbering.
     * @param grid_coordinates Cell location provided in terms of grid coordinates (i, j, k).
     * @throws std::out_of_range when the grid coordinates are outside of this grid subdivisions (nx, ny, nz).
     */
    virtual std::array<Index, NumberOfNodes> nodes(const VecInt & grid_coordinates) const;

    /** Get the position of the node node_id **/
    virtual VecFloat position(const Index & node_id) const;

    /**
     * Subdivide the cell at location cell_index
     * @param cell_index The index of the cell.
     * @throws std::logic_error when the cell to be subdivided is not a leaf-cell (is already subdivided).
     * @throws std::out_of_range when the cell_index doesn't point to any cell within this grid (subdivided cells indices included).
     */
    virtual void subdivide(const Index & cell_index);

    /** Get an regular hexahedron geometric representation (in space) of the given cell */
    geometry::RegularLinearHexahedron hexahedron_from_cell(const CellType & cell);

    /** Get an iterator over this grid leaf-cells. */
    CellRange leaf_cells() const {
        return CellRange(this);
    }

protected:
    ///< The anchor point position. It should be positioned at the center of the grid.
    VecFloat m_anchor;

    ///< Number of sub-cells in the x, y and z directions respectively.
    const VecInt m_number_of_subdivisions;

    ///<  Dimension of the grid from the anchor point in the x, y and z directions respectively.
    const VecFloat m_dimensions;

    ///< The cells this grid contains
    std::vector<CellType> m_cells;

    ///< Contain the number of cells in a given level of subdivision. If the grid contains only leaf-cells, than only
    ///< the level 0 of subdivision exists and the number cells of this level is the number of cells of this grid.
    ///< If two cells within this grid are subdivided, and one of them is again
    ///< subdivided, than the level 1 contains 2*CellType::NumberOfSubcells and the level 2 contains 1*CellType::NumberOfSubcells.
    ///< This counter will be used to assign numbers to every cell, nodes, edges and faces of this grid.
    std::vector<Int> m_number_of_cells_per_level;
};


template <class TCell>
struct Grid<TCell>::CellRange
{
    using GridType = Grid<TCell>;
    using CellType = TCell;
    using VecFloat = typename CellType::VecFloat;
    using Index = typename CellType::Index;
    using VecInt = typename CellType::VecInt;
    using Int = typename VecInt::ValueType;
    using Float = typename VecFloat::ValueType;

    struct CellIterator;

    explicit CellRange(const GridType * const grid) : p_grid(grid) {}

    CellIterator begin() const {
        return CellIterator(p_grid, 0, p_grid->get(0));
    }

private:
    const GridType * const p_grid;
};

template <class TCell>
struct Grid<TCell>::CellRange::CellIterator
{
    using GridType = Grid<TCell>;
    using CellType = TCell;
    using Int = typename CellType::VecInt::ValueType;

    explicit CellIterator(const GridType * const grid, const Int outer_cell_index, const CellType * current_cell)
    : p_grid(grid), p_outer_cell_index(outer_cell_index), p_current_cell(current_cell)
    {}

private:
    ///< The grid on which we are iterating over its leaf cells.
    const GridType * const p_grid;

    ///< This represent the index of the outer cell (in the level 0 of the grid) that contains p_current_cell
    ///< (which can be the outer cell itself if it is a leaf, or one of its sub-cells if it is subdivided)
    Int p_outer_cell_index;

    ///< The current cell (or sub-cell) pointing by the iterator. It will always be a leaf-cell.
    const CellType * p_current_cell;
};

} // namespace engine

} // namespace topology

} // namespace caribou

#endif //CARIBOU_TOPOLOGY_ENGINE_GRID_GRID_H
