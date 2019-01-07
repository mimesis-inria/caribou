#ifndef CARIBOU_TOPOLOGY_ENGINE_GRID_GRID_H
#define CARIBOU_TOPOLOGY_ENGINE_GRID_GRID_H

#include <vector>
#include <memory>
#include <cmath>

#include <Caribou/config.h>

#include <Caribou/Topology/Engine/Grid/Internal/UnidimensionalGrid.h>
#include <Caribou/Topology/Engine/Grid/Internal/MultidimensionalGrid.h>

namespace caribou {
namespace topology {
namespace engine {

template <size_t Dim, class CellType_=void>
struct Grid : public internal::BaseMultidimensionalGrid<Dim, Grid<Dim, CellType_>>
{
};

template <size_t Dim>
struct Grid<Dim, void> : public internal::BaseMultidimensionalGrid<Dim, Grid<Dim, void>>
{
    static constexpr size_t Dimension = Dim;

    using Base = internal::BaseMultidimensionalGrid<Dimension, Grid<Dimension, void>>;
    using Base::Base;
};

template <>
struct Grid<1, void> : public internal::BaseUnidimensionalGrid<Grid<1, void>>
{
    static constexpr size_t Dimension = 1;

    using Base = internal::BaseUnidimensionalGrid<Grid<Dimension, void>>;
    using Base::Base;
};

///**
// * A Grid is a set of multiple cell entities (same-length lines in 1D, rectangles
// * in 2D and rectangular hexahedrons in 3D) aligned in the x, y and z axis.
// */
//template <size_t Dim>
//struct Grid : public BaseGrid<Dim, Grid<Dim>>
//{
//    static constexpr size_t Dimension = Dim;
//
//    using Base = BaseGrid<Dim, Grid<Dim>>;
//
//    using Int = size_t;
//    using Float = FLOATING_POINT_TYPE;
//
//    using VecFloat = caribou::algebra::Vector<Dimension, Float>;
//    using VecInt = caribou::algebra::Vector<Dimension, Int>;
//
//    using Index = Int;
//    using NodeIndex = Int;
//    using CellIndex = Int;
//    using Position = VecFloat;
//    using Size = VecFloat;
//
//    /** Default constructor is not permitted **/
//    Grid() = delete;
//
//    /**
//     * Constructor of a grid
//     * @param anchor The anchor point position vector (x, y [, z]). It should be positioned at the center of the grid.
//     * @param subdivisions Vector of integers (nx, ny [, nz]) which specify the number of sub-cells in the x, y [and z]
//     *        directions respectively.
//     * @param dimensions Vector of floats (sx, sy [, sz]) which specify the size of the grid from the anchor point in
//     *        the x, y [and z] directions respectively.
//     */
//    Grid(const Position & anchor, const VecInt & subdivisions, const Size & dimensions)
//            : m_anchor(anchor), m_number_of_subdivisions(subdivisions), m_dimensions(dimensions)
//    {}
//
//    /** Get the number of cell subdivisions (nx, ny, nz) of this grid. **/
//    inline VecInt
//    number_of_subdivision() const
//    {return m_number_of_subdivisions;};
//
//    /** Get this grid dimensions (sx, sy, sz). **/
//    inline Size
//    size() const
//    { return m_dimensions;};
//
//    /** Get dimensions (hx, hy, hz) of a top-cell in this grid. **/
//    inline VecFloat
//    cell_size() const
//    {
//        const auto & N = number_of_subdivision();
//        const auto & S = size();
//        const auto   H = S.direct_division(N);
//
//        return H;
//    };
//
//    /** Get the position of the node node_id relative to the world frame (0, 0, 0). **/
//    inline Position
//    position(const NodeIndex & index) const
//    {
//        // Relative position of the node within the grid
//        const auto p = Base::relative_position(index);
//
//        return m_anchor + p; // World position
//    }
//
//protected:
//    ///< The anchor point position. It should be positioned at the center of the grid.
//    Position m_anchor;
//
//    ///< Number of sub-cells in the x, y and z directions respectively.
//    const VecInt m_number_of_subdivisions;
//
//    ///<  Dimension of the grid from the anchor point in the x, y and z directions respectively.
//    const Size m_dimensions;
//};
//
//
//template <class TCell>
//struct Grid<TCell>::CellRange
//{
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
//};
//
//template <class TCell>
//struct Grid<TCell>::CellRange::CellIterator
//{
//    using Self = Grid<TCell>::CellRange::CellIterator;
//    using GridType = Grid<TCell>;
//    using CellType = TCell;
//    using Int = typename CellType::VecInt::ValueType;
//
//    explicit CellIterator(const GridType * const grid, const Int outer_cell_index, const CellType * current_cell)
//    : p_grid(grid), p_outer_cell_index(outer_cell_index), p_current_cell(current_cell)
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

} // namespace engine
} // namespace topology
} // namespace caribou

#endif //CARIBOU_TOPOLOGY_ENGINE_GRID_GRID_H
