#ifndef CARIBOU_COMPONENTS_TOPOLOGY_MULTILEVELSPARSEGRIDTOPOLOGY_H
#define CARIBOU_COMPONENTS_TOPOLOGY_MULTILEVELSPARSEGRIDTOPOLOGY_H

#include <SofaBaseTopology/SparseGridTopology.h>
#include <array>

namespace caribou {
namespace components {
namespace topology {

using namespace sofa::core::objectmodel;

class MultilevelSparseGridTopology : public sofa::component::topology::SparseGridTopology
{
public:

    using Parent = sofa::component::topology::SparseGridTopology;
    using CellType = Parent::Type;
    using Vec3i = Parent::Vec3i;
    using HexaID = Parent::HexaID;

    struct Cell {
        unsigned char level; ///< If the level is 0, it is a cell, if the level is 1, it is the cell's subcell, if the level is 2, it is the subcell(1) subcell, etc
        CellType type; ///< Indicates if this cell is outside, inside or on the boundary.
        Cell * parent; ///< Parent cell (or subcell if parent->level > 0)
        Vec3i position; ///< The position of the cell (or subcell) in the parent context. If it is a cell, the position is relative to the grid. If it is a subcell, the position is relative to its parent cell.
        std::array<Cell*, 8> subcells; ///< List of subcells (there should be 8 per cell or subcell).

        Cell() : level(0), type(Type::INSIDE), parent(nullptr), position(), subcells() {}
        Cell (unsigned char level, CellType type, Vec3i position, Cell* parent)
                : level(level), type(type), parent(parent), position(position), subcells() {}
    };

    MultilevelSparseGridTopology();
    void init() override;

protected:
    Data<unsigned char> d_number_of_subdivision;

    sofa::helper::vector<Cell> m_cells;
};

} // namespace topology
} // namespace components
} // namespace caribou

#endif //CARIBOU_COMPONENTS_TOPOLOGY_MULTILEVELSPARSEGRIDTOPOLOGY_H
