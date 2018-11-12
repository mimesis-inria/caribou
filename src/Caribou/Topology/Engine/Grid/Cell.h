#ifndef CARIBOU_TOPOLOGY_ENGINE_GRID_CELL_H
#define CARIBOU_TOPOLOGY_ENGINE_GRID_CELL_H

#include <Caribou/config.h>
#include <Caribou/Algebra/Vector.h>

#include <memory>

namespace caribou
{

namespace topology
{

namespace engine
{

template <unsigned char Dim, class TCell>
struct Grid;

/**
 * A Cell represents a 2D quad (resp. 3D hexahedron) entity in a rectangular Grid made of multiple cells. It is always
 * contained in a parent grid, and it can also be subdivided in sub-cells by adding it a grid. Otherwise, if no grid is
 * added, the cell is said to be a leaf of its parent grid.
 *
 * @tparam Dimension The dimension (2D or 3D) of the cell.
 */
template <unsigned char Dim>
struct Cell
{
    static constexpr size_t Dimension = Dim;
    static_assert(Dimension == 2 or Dimension == 3, "A grid cell must be in two or three dimension.");

    using VecFloat = algebra::Vector<Dimension, FLOATING_POINT_TYPE>;
    using Index = size_t;
    using VecInt = algebra::Vector<Dimension, Index>;
    using GridType = Grid<Dimension, Cell<Dimension>>;

    static constexpr size_t NumberOfNodes = (unsigned char) (1 << Dimension);

    /** Default constructor is not permitted **/
    explicit Cell() = delete;

    Cell(GridType* parent, Index index) : m_grid(nullptr), m_parent(parent), m_index (index) {
        if (!parent) {
            throw std::logic_error("A cell must be contained in a parent grid.");
        }
    };

    /**
     * Subdivide the cell into nx, ny and nz sub-cells
     * @param subdivisions Specify the number (nx, ny, nz) of sub-cells
     * @throws std::logic_error When the cell is already subdivided.
     * @return The created grid containing the sub-cells
     */
    GridType *
    subdivide(VecInt subdivisions)
    {
        if (!is_a_leaf()) {
            throw std::logic_error("Trying to subdivide an already subdivided cell.");
        }

        Index anchor_index = nodes()[0];
        VecFloat anchor_position = m_parent->position(anchor_index);
        VecFloat cell_dimensions = size();

        m_grid.reset(new GridType(anchor_position, subdivisions, cell_dimensions));

        return m_grid.get();
    };

    /** This cell size (hx, hy, hz) **/
    inline VecFloat
    size() const {
        return m_parent->cell_size();
    };

    /** True if the cell is a leaf (it contains no sub-cells) **/
    inline bool
    is_a_leaf() const
    {return (m_grid == nullptr);};

    /** Get the index of the cell (relative to its parent grid). **/
    inline Index
    index() const {return m_index;};

    /** Get the indices of the nodes forming the cell. */
    inline std::array<Index, NumberOfNodes>
    nodes() const {
        return m_parent->nodes(m_parent->grid_coordinates(m_index));
    };

protected:
    std::unique_ptr<GridType> m_grid = nullptr;
    GridType* m_parent = nullptr;
    const Index m_index;
};

} // namespace engine

} // namespace topology

} // namespace caribou

#endif //CARIBOU_TOPOLOGY_ENGINE_GRID_CELL_H
