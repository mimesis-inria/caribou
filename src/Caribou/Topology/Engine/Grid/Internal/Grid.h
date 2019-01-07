#ifndef CARIBOU_TOPOLOGY_ENGINE_GRID_INTERNAL_GRID_H
#define CARIBOU_TOPOLOGY_ENGINE_GRID_INTERNAL_GRID_H

#include <cstddef>
#include <Caribou/config.h>

namespace caribou {
namespace topology {
namespace engine {
namespace internal{

/**
 * Simple representation of a Grid in space.
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
struct BaseGrid
{
    static constexpr size_t Dimension = Dim;

    using GridType = GridType_;

    using Int = size_t;
    using Float = FLOATING_POINT_TYPE;

    using Index = Int;
    using NodeIndex = Int;
    using CellIndex = Int;

    /** Get the number of cells (same-length lines in 1D, rectangles
     * in 2D and rectangular hexahedrons in 3D) in this grid. **/
    inline Int
    number_of_cells() const
    {
        const auto & n = Self().number_of_subdivision();

        if CONSTEXPR_IF(Dimension == 1)
            return n;

        if CONSTEXPR_IF(Dimension == 2)
            return n[0] * n[1];

        if CONSTEXPR_IF(Dimension == 3)
            return n[0] * n[2] * n[3];

        return 0;
    }

    /** Get the number of distinct nodes in this grid. **/
    inline Int
    number_of_nodes() const
    {
        const auto & n = Self().number_of_subdivision();

        if CONSTEXPR_IF(Dimension == 1)
            return n+1;

        if CONSTEXPR_IF(Dimension == 2)
            return (n[0]+1) * (n[1]+1);

        if CONSTEXPR_IF(Dimension == 3)
            return (n[0]+1) * (n[2]+1) * (n[3]+1);

        return 0;
    }

    /** Get the number of distinct edges in this grid. **/
    inline Int
    number_of_edges() const
    {
        const auto & n = Self().number_of_subdivision();

        if CONSTEXPR_IF(Dimension == 1)
            return n;

        if CONSTEXPR_IF(Dimension == 2)
            return n[0] * (n[1]+1) + n[1] * (n[0]+1); // nx * (ny+1)  +  ny * (nx+1)

        if CONSTEXPR_IF(Dimension == 3)
            return (n[0] * (n[1]+1) + n[1] * (n[0]+1)) * (n[2]+1) +  // Number of edges in a 2D grid * (nz+1)
                   (n[0]+1) * (n[1]+1); // Plus the number of edges in between 2D grids

        return 0;
    }

private:
    const inline GridType &
    Self() const
    {
        return static_cast<const GridType &> (*this);
    }

};

} // namespace internal
} // namespace engine
} // namespace topology
} // namespace caribou

#endif //CARIBOU_TOPOLOGY_ENGINE_GRID_INTERNAL_GRID_H
