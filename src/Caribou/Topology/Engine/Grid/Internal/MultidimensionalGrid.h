#ifndef CARIBOU_TOPOLOGY_ENGINE_GRID_INTERNAL_MULTIDIMENSIONALGRID_H
#define CARIBOU_TOPOLOGY_ENGINE_GRID_INTERNAL_MULTIDIMENSIONALGRID_H

#include <Caribou/Topology/Engine/Grid/Internal/Grid.h>
#include <Caribou/Algebra/Vector.h>

namespace caribou {
namespace topology {
namespace engine {
namespace internal {

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

    using Int = size_t;
    using Float = FLOATING_POINT_TYPE;

    using VecFloat = caribou::algebra::Vector<Dimension, Float>;
    using VecInt = caribou::algebra::Vector<Dimension, Int>;

    using Index = Int;
    using NodeIndex = Int;
    using CellIndex = Int;
    using RelativePosition = VecFloat;
    using Size = VecFloat;

    using BaseGrid<Dimension, GridType_>:: BaseGrid;

    /** Get the position of the node node_id relative to the grid frame (bottom-left most node). **/
    inline RelativePosition
    relative_position(const NodeIndex & index) const
    {
        const auto & n = Self().number_of_subdivision();
        const auto & nx = n[0]+1; // Number of nodes in the x direction
        const auto & ny = n[1]+1; // Number of nodes in the y direction

        const Index k = (Dimension == 3) ? (index / (nx*ny)) : 0; // Node indice in the z direction
        const Index j = (index - (k*nx*ny)) / nx; // Node indice in the y direction
        const Index i = index - ((k*nx*ny) + (j*nx)); // Node indice in the x direction

        const auto & h = Self().cell_size();
        const auto & hx = h[0]; // Width of a cell
        const auto & hy = h[1]; // Height of a cell
        const auto & hz = h[2]; // Dept of a cell

        // Relative position of the node within the grid
        const VecFloat p = (Dimension == 2) ? VecFloat({i*hx, j*hy}) : VecFloat({i*hx, j*hy, k*hz});

        return p; // Relative position
    }

    /** Query the cell that contains a given position */
    inline VecInt query(const VecFloat & position) const {

    }

private:
    inline constexpr
    const GridType &
    Self() const
    {
        return static_cast<const GridType &> (*this);
    }

};

} // namespace internal
} // namespace engine
} // namespace topology
} // namespace caribou

#endif //CARIBOU_TOPOLOGY_ENGINE_GRID_INTERNAL_MULTIDIMENSIONALGRID_H
