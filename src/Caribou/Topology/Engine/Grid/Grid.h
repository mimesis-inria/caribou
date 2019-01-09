#ifndef CARIBOU_TOPOLOGY_ENGINE_GRID_GRID_H
#define CARIBOU_TOPOLOGY_ENGINE_GRID_GRID_H

#include <vector>
#include <memory>
#include <cmath>

#include <Caribou/config.h>

#include <Caribou/Topology/Engine/Grid/Internal/UnidimensionalGrid.h>
#include <Caribou/Topology/Engine/Grid/Internal/MultidimensionalGrid.h>

namespace caribou {
namespace topology {
namespace engine {

template <size_t Dim, class CellType_=void>
struct Grid : public internal::BaseMultidimensionalGrid<Dim, Grid<Dim, CellType_>>
{
    using Base = internal::BaseMultidimensionalGrid<Dim, Grid<Dim, void>>;
    using Base::Base;
};

template <class CellType_>
struct Grid<1, CellType_> : public internal::BaseUnidimensionalGrid<Grid<1, CellType_>>
{
    static constexpr size_t Dimension = 1;

    using Base = internal::BaseUnidimensionalGrid<Grid<Dimension, void>>;
    using Base::Base;
};

} // namespace engine
} // namespace topology
} // namespace caribou

#endif //CARIBOU_TOPOLOGY_ENGINE_GRID_GRID_H
