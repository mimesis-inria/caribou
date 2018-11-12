#ifndef CARIBOU_TOPOLOGY_ENGINE_GRID_GRID_INL
#define CARIBOU_TOPOLOGY_ENGINE_GRID_GRID_INL

#include <Caribou/Topology/Engine/Grid/Grid.h>

namespace caribou
{

using namespace geometry;

namespace topology
{

namespace engine
{
// ########################
// Cell template functions
// ########################
template <unsigned char Dimension>
Grid<Dimension>*
Cell<Dimension>::subdivide (VecInt subdivisions)
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

template <unsigned char Dimension>
typename Cell<Dimension>::VecFloat
Cell<Dimension>::size () const
{
    return m_parent->cell_size();
};

template <unsigned char Dimension>
std::array<typename Cell<Dimension>::Index, Cell<Dimension>::NumberOfNodes>
Cell<Dimension>::nodes () const
{
    return m_parent->nodes(m_parent->grid_coordinates(m_index));
};

// ########################
// Grid template functions
// ########################
template <unsigned char Dimension>
typename Grid<Dimension>::CellType &
Grid<Dimension>::get (const VecInt &grid_coordinates)
{
    return *(cells[cell_index(grid_coordinates)]);
};

} // namespace engine

} // namespace topology

} // namespace caribou

#endif //CARIBOU_TOPOLOGY_ENGINE_GRID_GRID_INL