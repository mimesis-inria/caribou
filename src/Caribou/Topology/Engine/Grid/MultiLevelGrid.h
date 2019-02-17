#ifndef CARIBOU_TOPOLOGY_ENGINE_GRID_MULTILEVELGRID_H
#define CARIBOU_TOPOLOGY_ENGINE_GRID_MULTILEVELGRID_H

#include <vector>
#include <memory>

#include <Caribou/Topology/Engine/Grid/Grid.h>
#include <Caribou/config.h>

namespace caribou::topology::engine {

/**
 * A Multilevel grid is Grid (see caribou::topology::engine::Grid) that allows its cell to be split in a set of orthants
 * (2 segments in 1D, 4 quadrants in 2D and 8 octants in 3D).
 * @tparam TCell The type of cell this grid contains. It must inherits caribou::topology::engine::Cell
 */
template <size_t Dim, class TCell>
struct MultiLevelGrid : public Grid<Dim>
{
    using Base = Grid<Dim>;

    using CellType = TCell;

    static constexpr size_t Dimension = Dim;

    using Base::Base;
    using NodeIndex = typename Base::NodeIndex;
    using CellIndex = typename Base::CellIndex;
    using Dimensions = typename Base::Dimensions;
    using Subdivisions = typename Base::Subdivisions;
    using LocalCoordinates = typename Base::LocalCoordinates;
    using WorldCoordinates = typename Base::WorldCoordinates;
    using GridCoordinates = typename Base::GridCoordinates;
    using CellSet = typename Base::CellSet;
    using Base::Base;

    struct CellRange;

    /**
     * Get the cell located at grid_index (i, j, k).
     * @param grid_coordinates Cell location provided in terms of grid coordinates (i, j, k).
     */
    inline CellType &
    get(const GridCoordinates & coordinates) noexcept
    {
        return get(cell_index_at(coordinates));
    }

    /**
     * Get the cell located at grid_index (i, j, k).
     * @param grid_coordinates Cell location provided in terms of grid coordinates (i, j, k).
     */
    inline const CellType &
    get(const GridCoordinates & coordinates) const noexcept
    {
        return get(cell_index_at(coordinates));
    }

    /**
     * Get the cell located at a given index.
     */
    inline CellType &
    get(const CellIndex & index) noexcept
    {
        return m_cells[index];
    }

    /**
     * Get the cell located at a given index.
     */
    inline const CellType &
    get(const CellIndex & index) const noexcept
    {
        return m_cells[index];
    }

    /** Get the cell directly at the right (x axis) of the given cell. If there are no cells at the right (we reach the
     * boundary of the grid), null is returned. If there is multiple cells, only the cell of the same level (therefore
     * of the same dimension) as the one given is returned. If there are no cell at the same level of the given one, the
     * cell at the closest level is returned (this cell will therefore be larger than the given one). */
//    CellType * right_cell_of(const CellType & cell);

    /** Get an iterator over this grid leaf-cells. */
    CellRange leaf_cells() const {
//        return CellRange(this);
    }

protected:

    ///< The cells this grid contains
    std::vector<CellType> m_cells;

    ///< Contain the number of cells in a given level of subdivision. If the grid contains only leaf-cells, than only
    ///< the level 0 of subdivision exists and the number cells of this level is the number of cells of this grid.
    ///< If two cells within this grid are subdivided, and one of them is again
    ///< subdivided, than the level 1 contains 2*CellType::NumberOfSubcells and the level 2 contains 1*CellType::NumberOfSubcells.
    ///< This counter will be used to assign numbers to every cell, nodes, edges and faces of this grid.
//    std::vector<Int> m_number_of_cells_per_level;
};


template <class TCell>
struct Grid<TCell>::CellRange
{
//    using GridType = Grid<TCell>;
//    using CellType = TCell;
//    using VecFloat = typename CellType::VecFloat;
//    using Index = typename CellType::Index;
//    using VecInt = typename CellType::VecInt;
//    using Int = typename VecInt::ValueType;
//    using Float = typename VecFloat::ValueType;
//
//    struct CellIterator;
//
//    explicit CellRange(const GridType * const grid) : p_grid(grid) {}
//
//    CellIterator begin() const {
//        // Get the bottom-left most outer-cell of the grid
//        const CellType * bottom_left_cell = &p_grid->get(0);
//
//        // Go down its sub-cells tree until we get its bottom-left most leaf sub-cell
//        while (not bottom_left_cell->is_a_leaf()) {
//            bottom_left_cell = &bottom_left_cell->child(0);
//        }
//
//        return CellIterator(p_grid, 0, bottom_left_cell);
//    }
//
//    CellIterator end() const {
//        return CellIterator(p_grid, CellType::Nx * CellType::Ny * CellType::Nz, nullptr);
//    }
//
//private:
//    const GridType * const p_grid;
};

//template <class TCell>
//struct Grid<TCell>::CellRange::CellIterator
//{
//    using Self = Grid<TCell>::CellRange::CellIterator;
//    using GridType = Grid<TCell>;
//    using CellType = TCell;
//    using Int = typename CellType::VecInt::ValueType;
//
//    explicit CellIterator(const GridType * const grid, const Int outer_cell_index, const CellType * current_cell)
//            : p_grid(grid), p_outer_cell_index(outer_cell_index), p_current_cell(current_cell)
//    {}
//
//    Self& operator++() {
////        num = TO >= FROM ? num + 1: num - 1;
////        return *this;
//    }
//
//    Self operator++(int) {
//        Self retval = *this;
//        ++(*this);
//        return retval;
//    }
//
//    bool operator==(Self /*other*/) const {
////        return num == other.num;
//        return true;
//    }
//
//    bool operator!=(Self other) const {
//        return !(*this == other);
//    }
//
//    const CellType * operator*() const {
//        return p_current_cell;
//    }
//
//private:
//    ///< The grid on which we are iterating over its leaf cells.
//    const GridType * const p_grid;
//
//    ///< This represent the index of the outer cell (in the level 0 of the grid) that contains p_current_cell
//    ///< (which can be the outer cell itself if it is a leaf, or one of its sub-cells if it is subdivided)
//    Int p_outer_cell_index;
//
//    ///< The current cell (or sub-cell) pointing by the iterator. It will always be a leaf-cell.
//    const CellType * p_current_cell;
//};

} // namespace caribou::topology::engine

#endif //CARIBOU_TOPOLOGY_ENGINE_GRID_MULTILEVELGRID_H
