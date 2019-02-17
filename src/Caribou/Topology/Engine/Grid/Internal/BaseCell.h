#ifndef CARIBOU_TOPOLOGY_ENGINE_GRID_INTERNAL_BASECELL_H
#define CARIBOU_TOPOLOGY_ENGINE_GRID_INTERNAL_BASECELL_H

#include <Caribou/config.h>
#include <Caribou/Algebra/Vector.h>

#include <memory>

namespace caribou::topology::engine::internal {

/**
 * Simple representation of a Cell.
 *
 * ** Do not use this class directly. Use instead caribou::topology::engine::Cell. **
 *
 * The functions declared in this class can be used with any type of cells.
 *
 * Do to so, it uses the Curiously recurring template pattern (CRTP) :
 *    https://en.wikipedia.org/wiki/Curiously_recurring_template_pattern
 *
 * A cell represents a segment in 1D, a rectangle quad in 2D or a 3D rectangle hexahedron in 3D.
 *
 * @tparam Dim The dimension of the cell (1D, 2D or 3D).
 * @tparam CellType_ Type of the derived cell class that will implement the final functions.
 */
template <unsigned char Dim, class CellType_>
struct BaseCell
{
    static constexpr char Dimension = Dim;
    static_assert(Dimension == 1 or Dimension == 2 or Dimension == 3, "A cell must be in one, two or three dimension.");

    BaseCell() = delete;

    template <typename GridType, typename IndexType>
    BaseCell(const GridType & /*grid*/, const IndexType & /*cell_index*/) {}
};

} // namespace caribou::topology::engine::internal

#endif //CARIBOU_TOPOLOGY_ENGINE_GRID_INTERNAL_BASECELL_H
