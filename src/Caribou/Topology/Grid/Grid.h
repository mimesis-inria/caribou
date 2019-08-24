#ifndef CARIBOU_TOPOLOGY_ENGINE_GRID_GRID_H
#define CARIBOU_TOPOLOGY_ENGINE_GRID_GRID_H

#include <vector>
#include <memory>
#include <cmath>
#include <array>

#include <Caribou/config.h>
#include <Caribou/Geometry/Traits.h>

#include <Caribou/Topology/Grid/Internal/BaseUnidimensionalGrid.h>
#include <Caribou/Topology/Grid/Internal/BaseMultidimensionalGrid.h>

#include <Caribou/Geometry/Segment.h>
#include <Caribou/Geometry/Quad.h>
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
    using Element = geometry::Segment<1, geometry::interpolation::Segment2>;

    inline
    Element
    cell_at(const GridCoordinates & coordinates) const
    {
        return cell_at(cell_index_at(coordinates));
    }

    inline
    Element
    cell_at(const CellIndex & index) const {
        const auto n0 = node(index);
        const auto n1 = node(index + 1);
        return Element (n0, n1);
    }

    inline
    std::array<NodeIndex, caribou::traits<Element>::NumberOfNodes>
    node_indices_of(const GridCoordinates & coordinates) const
    {
        return node_indices_of(cell_index_at(coordinates));
    }

    inline
    std::array<NodeIndex, caribou::traits<Element>::NumberOfNodes>
    node_indices_of(const CellIndex & index) const
    {
        return {{
            index,
            index + 1
        }};
    }
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
    using Element = geometry::Quad<2, geometry::interpolation::Quad4>;

    inline
    Element
    cell_at(const GridCoordinates & cell_coordinates) const
    {
        GridCoordinates e1; e1 << 1, 0;
        GridCoordinates e2; e2 << 0, 1;

        const auto n0 = node(node_index_at(cell_coordinates));
        const auto n1 = node(node_index_at(cell_coordinates + e1));
        const auto n2 = node(node_index_at(cell_coordinates + e1 + e2));
        const auto n3 = node(node_index_at(cell_coordinates + e2));

        return Element(n0, n1, n2, n3);
    }

    inline
    Element
    cell_at(const CellIndex & cell_index) const
    {
        return cell_at(grid_coordinates_at(cell_index));
    }

    inline
    std::array<NodeIndex, caribou::traits<Element>::NumberOfNodes>
    node_indices_of(const GridCoordinates & cell_coordinates) const
    {
        GridCoordinates e1; e1 << 1, 0;
        GridCoordinates e2; e2 << 0, 1;
        
        return {{
            node_index_at(cell_coordinates), 
            node_index_at(cell_coordinates + e1),
            node_index_at(cell_coordinates + e1 + e2),
            node_index_at(cell_coordinates + e2)
        }};
    }

    inline
    std::array<NodeIndex, caribou::traits<Element>::NumberOfNodes>
    node_indices_of(const CellIndex & index) const
    {
        return node_indices_of(grid_coordinates_at(index));
    }
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
    using Element = geometry::RectangularHexahedron<geometry::interpolation::Hexahedron8>;

    inline
    Element
    cell_at(const GridCoordinates & coordinates) const
    {
        const auto H = Base::H();
        const auto center = Base::m_anchor_position + (coordinates.array().cast<FLOATING_POINT_TYPE>() * H.array()).matrix() + H/2.;
        return Element (center, H);
    }

    inline
    Element
    cell_at(const CellIndex & cell_index) const
    {
        return cell_at(grid_coordinates_at(cell_index));
    }

    inline
    std::array<NodeIndex, caribou::traits<Element>::NumberOfNodes>
    node_indices_of(const GridCoordinates & cell_coordinates) const
    {
        GridCoordinates e1; e1 << 1, 0, 0;
        GridCoordinates e2; e2 << 0, 1, 0;
        GridCoordinates e3; e3 << 0, 0, 1;

        return {{
            node_index_at(cell_coordinates),
            node_index_at(cell_coordinates + e1),
            node_index_at(cell_coordinates + e1 + e2),
            node_index_at(cell_coordinates + e2),
            node_index_at(cell_coordinates + e3),
            node_index_at(cell_coordinates + e1 + e3),
            node_index_at(cell_coordinates + e1 + e2 + e3),
            node_index_at(cell_coordinates + e2 + e3)
        }};
    }

    inline
    std::array<NodeIndex, caribou::traits<Element>::NumberOfNodes>
    node_indices_of(const CellIndex & index) const
    {
        return node_indices_of(grid_coordinates_at(index));
    }

    /** Get the number of distinct edges in this grid. **/
    inline UInt
    number_of_faces() const noexcept
    {
        const auto & N = Base::N();

        const auto & nx = N[0];
        const auto & ny = N[1];
        const auto & nz = N[2];

        const auto number_of_faces_in_2D_grid = nx * ny;
        const auto number_of_faces_between_2D_grids = (nx * ny) + ((nx + 1) * ny) + (nx * (ny + 1));

        return number_of_faces_in_2D_grid * (nz+1)
               + number_of_faces_between_2D_grids * nz;
    }
};

} // namespace caribou::topology

#endif //CARIBOU_TOPOLOGY_ENGINE_GRID_GRID_H
