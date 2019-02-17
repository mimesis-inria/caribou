#ifndef CARIBOU_TOPOLOGY_ENGINE_GRID_GRIDCONTAINER_H
#define CARIBOU_TOPOLOGY_ENGINE_GRID_GRIDCONTAINER_H

#include <vector>
#include <memory>

#include <Caribou/Topology/Engine/Grid/Grid.h>
#include <Caribou/config.h>

namespace caribou::topology::engine {

/**
 * A grid container is Grid (see caribou::topology::engine::Grid) that keeps in memory a given structure for each
 * of its cells. The cell type is defined by the template argument TCell.
 *
 * @tparam TCell The type of cell this grid contains.
 */
template<size_t Dim, class TCell>
struct GridContainer : public Grid<Dim>
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

    struct CellRange;

    /**
     * Get the cell located at grid_index (i, j, k).
     * @param grid_coordinates Cell location provided in terms of grid coordinates (i, j, k).
     */
    inline CellType &
    get(const GridCoordinates &coordinates) noexcept
    {
        return get(cell_index_at(coordinates));
    }

    /**
     * Get the cell located at grid_index (i, j, k).
     * @param grid_coordinates Cell location provided in terms of grid coordinates (i, j, k).
     */
    inline const CellType &
    get(const GridCoordinates &coordinates) const noexcept
    {
        return get(cell_index_at(coordinates));
    }

    /**
     * Get the cell located at a given index.
     */
    inline CellType &
    get(const CellIndex &index) noexcept
    {
        return m_cells[index];
    }

    /**
     * Get the cell located at a given index.
     */
    inline const CellType &
    get(const CellIndex &index) const noexcept
    {
        return m_cells[index];
    }

protected:

    ///< The cells this grid contains
    std::vector<CellType> m_cells;
};

} // namespace caribou::topology::engine
#endif //CARIBOU_TOPOLOGY_ENGINE_GRID_GRIDCONTAINER_H
