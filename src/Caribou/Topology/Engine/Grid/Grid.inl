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

template <char Dimension>
Grid<Dimension>*
Cell<Dimension>::subdivide (VecInt subdivisions)
{
    if (!is_a_leaf()) {
        throw std::logic_error("Trying to subdivide an already subdivided cell.");
    }

};

template <char Dimension>
typename Cell<Dimension>::VecFloat
Cell<Dimension>::size () const
{
    return m_parent->cell_size();
};

template <char Dimension>
Cell<Dimension> &
Grid<Dimension>::get (const VecInt &grid_coordinates)
{
    return *(cells[cell_index(grid_coordinates)]);
};

} // namespace engine

} // namespace topology

} // namespace caribou

#endif //CARIBOU_TOPOLOGY_ENGINE_GRID_GRID_INL