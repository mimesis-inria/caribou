#ifndef CARIBOU_TOPOLOGY_ENGINE_GRID_INTERNAL_UNIDIMENSIONALGRID_H
#define CARIBOU_TOPOLOGY_ENGINE_GRID_INTERNAL_UNIDIMENSIONALGRID_H

#include <Caribou/Topology/Engine/Grid/Internal/Grid.h>

namespace caribou {
namespace topology {
namespace engine {
namespace internal {

/**
 * Simple representation of an unidimensional (1D) Grid in space.
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
 * @tparam GridType_ Type of the derived grid class that will implement the final functions.
 */
template <class GridType_>
struct BaseUnidimensionalGrid : public BaseGrid<1, GridType_>
{
    static constexpr size_t Dimension = 1;

    using GridType = GridType_;

    using Int = size_t;
    using Float = FLOATING_POINT_TYPE;

    using Index = Int;
    using NodeIndex = Int;
    using CellIndex = Int;
    using RelativePosition = Float;
    using Size = Float;

    /** Default constructor **/
    constexpr
    BaseUnidimensionalGrid() = delete;

    /**
     * Constructor of a 1D grid (which is a consecutive set of same-length 1D segments along a line).
     * @param subdivisions Number of sub-cells (same-length 1D segments) in the direction of the line.
     * @param dimension Size of grid line from the anchor.
     */
    constexpr
    BaseUnidimensionalGrid(const Int & subdivisions, const Size & dimensions)
            : m_number_of_subdivisions(subdivisions), m_dimensions(dimensions)
    {}

    inline constexpr
    Float
    cell_size() const
    {
        return m_dimensions / m_number_of_subdivisions;
    }

    /** Get the position of the node node_id relative to the grid frame (bottom-left most node). **/
    inline constexpr
    RelativePosition
    relative_position(const NodeIndex & index) const
    {

        return index * Self().cell_size(); // Relative position
    }

protected:
    ///< Number of sub-cells in the x, y and z directions respectively.
    const Int m_number_of_subdivisions;

    ///<  Dimension of the grid from the anchor point in the x, y and z directions respectively.
    const Size m_dimensions;

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

#endif //CARIBOU_TOPOLOGY_ENGINE_GRID_INTERNAL_UNIDIMENSIONALGRID_H
