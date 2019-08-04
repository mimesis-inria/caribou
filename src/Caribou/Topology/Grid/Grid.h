#ifndef CARIBOU_TOPOLOGY_ENGINE_GRID_GRID_H
#define CARIBOU_TOPOLOGY_ENGINE_GRID_GRID_H

#include <vector>
#include <memory>
#include <cmath>

#include <Caribou/config.h>

#include <Caribou/Topology/Grid/Internal/BaseUnidimensionalGrid.h>
#include <Caribou/Topology/Grid/Internal/BaseMultidimensionalGrid.h>

#include <Caribou/Geometry/RectangularHexahedron.h>

namespace caribou::topology {

template <size_t Dim>
struct Grid
{
    static_assert(Dim > 0 and Dim < 4, "A grid can only be of dimension 1, 2 or 3.");
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

template <>
struct Grid<2> : public internal::BaseMultidimensionalGrid<2, Grid<2>>
{
    static constexpr size_t Dimension = 2;

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
struct Grid<3> : public internal::BaseMultidimensionalGrid<3, Grid<3>>
{
    static constexpr size_t Dimension = 3;

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

    inline
    geometry::RectangularHexahedron<geometry::interpolation::Hexahedron8>
    hexahedron(const GridCoordinates & coordinates) const
    {
        const auto H = Base::H();
        const auto center = Base::m_anchor_position + (coordinates.array().cast<FLOATING_POINT_TYPE>() * H.array()).matrix() + H/2.;
        return geometry::RectangularHexahedron<geometry::interpolation::Hexahedron8> (center, H);
    }
};

} // namespace caribou::topology

#endif //CARIBOU_TOPOLOGY_ENGINE_GRID_GRID_H
