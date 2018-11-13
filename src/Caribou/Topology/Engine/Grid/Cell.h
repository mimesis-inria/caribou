#ifndef CARIBOU_TOPOLOGY_ENGINE_GRID_CELL_H
#define CARIBOU_TOPOLOGY_ENGINE_GRID_CELL_H

#include <Caribou/config.h>
#include <Caribou/Algebra/Vector.h>
#include <Caribou/Geometry/Hexahedron.h>

#include <memory>

namespace caribou
{

using namespace geometry;

namespace topology
{

namespace engine
{

/**
 * A Cell represents a 2D rectangle quad (resp. 3D rectangle hexahedron) entity. It can be
 * contained in a parent cell, and it can also be subdivided in 4 sub-cells in 2D (resp. 8 sub-cells in 3D).
 * Otherwise, if no subdivision is done, the cell is said to be a leaf cell.
 *
 * @tparam Dimension The dimension (2D or 3D) of the cell.
 */
template <unsigned char Dim>
struct Cell
{
    static constexpr size_t Dimension = Dim;
    static_assert(Dimension == 2 or Dimension == 3, "A cell must be in two or three dimension.");

    static constexpr size_t NumberOfNodes = (unsigned char) (1 << Dimension);
    static constexpr size_t NumberOfSubcells = (unsigned char) (1 << Dimension);

    using VecFloat = algebra::Vector<Dimension, FLOATING_POINT_TYPE>;
    using Index = size_t;
    using VecInt = algebra::Vector<Dimension, Index>;

    using Subcells = std::array<Cell<Dimension>, NumberOfSubcells>;

    /** Default constructor */
    Cell() {};

    /**
     * Constructor with index assignment
     * @param index The index of this cell in its parent container
     */
    Cell(Index index) : m_index (index) {};

    /**
     * Subdivide the cell into 4 sub-cells in 2D (resp. 8 sub-cells in 3D)
     * @throws std::logic_error When the cell is already subdivided.
     * @return The current cell containing the sub-cells
     */
    Cell<Dimension> &
    subdivide()
    {
        if (!is_a_leaf()) {
            throw std::logic_error("Trying to subdivide an already subdivided cell.");
        }

        m_subcells.reset(new Subcells);
        for (Index subcell_index = 0; subcell_index < NumberOfSubcells; ++subcell_index) {
            Cell<Dimension> & cell = (*m_subcells)[subcell_index];
            cell.m_parent = this;
            cell.m_index = subcell_index;
            cell.m_level = m_level+1;
        }

        return *this;
    };

    /** True if the cell is a leaf (it contains no sub-cells) **/
    inline bool
    is_a_leaf() const
    {return (m_subcells == nullptr);};

    /** True if the cell got a parent cell **/
    inline bool
    has_parent() const
    {return (m_parent != nullptr);};

    /** Get the index of the cell (relative to its parent cell). **/
    inline Index
    index() const {return m_index;};

    /** Get the index of the cell (relative to its parent cell). **/
    inline Index
    level() const {return m_level;};

    /** Get the parent cell. **/
    inline const Cell<Dimension> &
    parent() const
    {
        if (not m_parent)
            throw std::logic_error("This cell does not contain a parent cell.");
        return *m_parent;
    };

    /** Get the child cell at specified index. **/
    inline Cell<Dimension> &
    child(Index index)
    {
        if (is_a_leaf())
            throw std::logic_error("This cell is a leaf cell, hence it does not contain any child cells.");

        if (index >= NumberOfSubcells)
            throw std::logic_error(
                    "Trying to access the " + std::to_string(index + 1)+ "th child cell, but this cell only have "
                    + std::to_string(NumberOfSubcells) + " child cells."
            );

        return *((*m_subcells)[index]);
    };

    /** Get the child cell at specified index. **/
    inline const Cell<Dimension> &
    child(Index index) const
    { return child(index); };

    /** Get the child cell at specified grid coordinate. **/
    inline Cell<Dimension> &
    child (const VecInt & grid_coordinates)
    {
        constexpr size_t nx = 2;

        const typename VecInt::ValueType i = grid_coordinates[0];
        const typename VecInt::ValueType j = grid_coordinates[1];

        // Note: The following condition block will be simplified by the compiler
        // since Dimension is a constexpr

        if (Dimension == 2) {
            // Dimension == 2
            return child(j*nx + i);
        } else {
            // Dimension == 3
            constexpr size_t ny = 2;
            const typename VecInt::ValueType k = grid_coordinates[2];
            return child(k*ny*nx + j*nx + i);
        }
    };

    /** Get the child cell at specified grid coordinate. **/
    inline const Cell<Dimension> &
    child (const VecInt & grid_coordinates) const
    { return child(grid_coordinates); };

protected:
    std::unique_ptr<Subcells> m_subcells = nullptr;
    Cell<Dimension> * m_parent = nullptr;
    Index m_index = 0;
    Index m_level = 0;
};

} // namespace engine

} // namespace topology

} // namespace caribou

#endif //CARIBOU_TOPOLOGY_ENGINE_GRID_CELL_H
