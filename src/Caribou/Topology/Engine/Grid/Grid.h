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

template <size_t Dim>
struct Grid : public internal::BaseMultidimensionalGrid<Dim, Grid<Dim>>
{
    static constexpr size_t Dimension = Dim;

    using Base = internal::BaseMultidimensionalGrid<Dimension, Grid<Dimension>>;
    using Base::Base;

    using NodeIndex = typename Base::NodeIndex;
    using CellIndex = typename Base::CellIndex;
    using Dimensions = typename Base::Dimensions;
    using Subdivisions = typename Base::Subdivisions;
    using LocalCoordinates = typename Base::LocalCoordinates;
    using WorldCoordinates = typename Base::WorldCoordinates;
    using GridCoordinates = typename Base::GridCoordinates;
    using CellSet = std::list<CellIndex>;
};

template <>
struct Grid<1> : public internal::BaseUnidimensionalGrid<Grid<1>>
{
    static constexpr size_t Dimension = 1;

    using Base = internal::BaseUnidimensionalGrid<Grid<Dimension>>;
    using Base::Base;

    using NodeIndex = typename Base::NodeIndex;
    using CellIndex = typename Base::CellIndex;
    using Dimensions = typename Base::Dimensions;
    using Subdivisions = typename Base::Subdivisions;
    using LocalCoordinates = typename Base::LocalCoordinates;
    using WorldCoordinates = typename Base::WorldCoordinates;
    using GridCoordinates = typename Base::GridCoordinates;
    using CellSet = std::list<CellIndex>;
};

} // namespace engine
} // namespace topology
} // namespace caribou

#endif //CARIBOU_TOPOLOGY_ENGINE_GRID_GRID_H
