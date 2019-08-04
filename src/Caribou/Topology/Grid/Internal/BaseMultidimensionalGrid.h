#ifndef CARIBOU_TOPOLOGY_GRID_INTERNAL_MULTIDIMENSIONALGRID_H
#define CARIBOU_TOPOLOGY_GRID_INTERNAL_MULTIDIMENSIONALGRID_H

#include <Caribou/Topology/Grid/Internal/BaseGrid.h>

namespace caribou::topology::internal {

/**
 * Simple representation of a multidimensional (2D or 3D) Grid in space.
 *
 * ** Do not use this class directly. Use instead caribou::topology::engine::Grid. **
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
struct BaseMultidimensionalGrid : public BaseGrid<Dim, GridType_>
{
    static constexpr size_t Dimension = Dim;

    static_assert(Dimension == 2 or Dimension == 3, "Only grid of dimension 2 or 3 is allowed.");

    using GridType = GridType_;
    using Base = BaseGrid<Dim, GridType_>;

    using NodeIndex = typename Base::NodeIndex;
    using CellIndex = typename Base::CellIndex;
    using Dimensions = typename Base::Dimensions;
    using Subdivisions = typename Base::Subdivisions;
    using LocalCoordinates = typename Base::LocalCoordinates;
    using WorldCoordinates = typename Base::WorldCoordinates;
    using GridCoordinates = typename Base::GridCoordinates;
    using CellSet = typename Base::CellSet;

    using Base::Base;

private:
    inline constexpr
    const GridType &
    Self() const
    {
        return static_cast<const GridType &> (*this);
    }

};

} // namespace caribou::topology::internal

#endif //CARIBOU_TOPOLOGY_GRID_INTERNAL_MULTIDIMENSIONALGRID_H
