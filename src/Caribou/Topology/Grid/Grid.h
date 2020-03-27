#ifndef CARIBOU_TOPOLOGY_GRID_GRID_H
#define CARIBOU_TOPOLOGY_GRID_GRID_H

#include <vector>
#include <memory>
#include <cmath>
#include <array>

#include <Caribou/config.h>
#include <Caribou/Geometry/Traits.h>

#include <Caribou/Topology/Grid/Internal/BaseUnidimensionalGrid.h>
#include <Caribou/Topology/Grid/Internal/BaseMultidimensionalGrid.h>

#include <Caribou/Geometry/Segment.h>
#include <Caribou/Geometry/RectangularQuad.h>
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

    using Element = geometry::Segment<1, geometry::interpolation::Segment2>;
    using ElementIndices = std::array<NodeIndex, Element::NumberOfNodes>;

    [[nodiscard]] inline
    auto
    cell_at(const GridCoordinates & coordinates) const -> Element
    {
        return cell_at(cell_index_at(coordinates));
    }

    [[nodiscard]] inline
    auto
    cell_at(const CellIndex & index) const -> Element {
        const auto n0 = node(index);
        const auto n1 = node(index + 1);
        return Element (n0, n1);
    }

    [[nodiscard]] inline
    auto
    node_indices_of(const GridCoordinates & coordinates) const -> ElementIndices
    {
        return node_indices_of(cell_index_at(coordinates));
    }

    [[nodiscard]] inline
    auto
    node_indices_of(const CellIndex & index) const -> ElementIndices
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

    using Element = geometry::RectangularQuad<2, geometry::interpolation::Quad4>;
    using Edge = Element::BoundaryType;

    using ElementNodes = std::array<NodeIndex, Element::NumberOfNodes>;
    using EdgeNodes = std::array<NodeIndex, Edge::NumberOfNodes>;

    [[nodiscard]] inline auto
    cell_at(const GridCoordinates & coordinates) const -> Element
    {
        const auto H = Base::H();
        const auto center = Base::m_anchor_position + (coordinates.array().cast<FLOATING_POINT_TYPE>() * H.array()).matrix() + H/2.;
        return Element (center, H);
    }

    [[nodiscard]] inline auto
    cell_at(const CellIndex & cell_index) const -> Element
    {
        return cell_at(cell_coordinates_at(cell_index));
    }

    [[nodiscard]] inline auto
    node_indices_of(const GridCoordinates & cell_coordinates) const -> ElementNodes
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

    [[nodiscard]] inline auto
    node_indices_of(const CellIndex & index) const -> ElementNodes
    {
        return node_indices_of(cell_coordinates_at(index));
    }

    /** Get the number of distinct edges in this grid. **/
    [[nodiscard]] inline auto
    number_of_edges() const noexcept -> UInt
    {
        // If the grid is
        // Y
        // `
        // `_ _ _ _
        // |_|_|_|_|
        // |_|_|_|_|
        // |_|_|_|_|... X
        // Then one row is
        // |_|_|_|_|
        // And the total number of edges is
        // 3 * |_|_|_|_|  +  _ _ _ _

        const auto & n = m_number_of_subdivisions;

        const auto & nx = n[0];
        const auto & ny = n[1];

        const auto nb_edges_per_row = (2*nx + 1);
        const auto nb_edges_total = (nb_edges_per_row * ny) + nx;

        return nb_edges_total;
    }

    /** Get the node indices of the edge at edge index */
    [[nodiscard]] inline auto
    edge(const EdgeIndex & index) const noexcept -> EdgeNodes
    {
        const auto & n = m_number_of_subdivisions;

        const auto & nx = n[0];

        // If the grid is
        // Y
        // `
        // `_ _ _ _
        // |_|_|_|_|
        // |_|_|_|_|
        // |_|_|_|_|... X
        // Then one row is
        // |_|_|_|_|
        // And the total number of edges is
        // 3 * |_|_|_|_|  +  _ _ _ _

        const auto nb_edges_per_row = (2*nx + 1);

        NodeIndex n0, n1;

        const auto row = index / nb_edges_per_row;
        const auto index_on_row = index % nb_edges_per_row;

        // If the row is
        // |_|_|_|_|
        // Then the edge is
        //  _ _ _ _  => x direction
        // | | | | | => y direction
        const bool edge_is_on_the_x_direction = index_on_row < nx;

        if (edge_is_on_the_x_direction) {
            n0 = row * (nx+1) + index_on_row;
            n1 = n0 + 1;
        } else {
            n0 = row * (nx+1) + index_on_row - nx;
            n1 = n0 + (nx + 1);
        }

        return {n0, n1};
    }
};

template <>
struct Grid<3> : public internal::BaseMultidimensionalGrid<3, Grid<3>>
{
    static constexpr size_t Dimension = 3;

    using Base = internal::BaseMultidimensionalGrid<Dimension, Grid<Dimension>>;
    using Base::Base;

    using FaceIndex = typename Base::UInt;

    using Element = geometry::RectangularHexahedron<geometry::interpolation::Hexahedron8>;
    using Face = Element::BoundaryType;
    using Edge = Face::BoundaryType;

    using ElementNodes = std::array<NodeIndex, Element::NumberOfNodes>;
    using FaceNodes = std::array<NodeIndex, Face::NumberOfNodes>;
    using EdgeNodes = std::array<NodeIndex, Edge::NumberOfNodes>;


    [[nodiscard]] inline auto
    cell_at(const GridCoordinates & coordinates) const noexcept -> Element
    {
        const auto H = Base::H();
        const auto center = Base::m_anchor_position + (coordinates.array().cast<FLOATING_POINT_TYPE>() * H.array()).matrix() + H/2.;
        return Element (center, H);
    }

    [[nodiscard]] inline auto
    cell_at(const CellIndex & cell_index) const noexcept -> Element
    {
        return cell_at(cell_coordinates_at(cell_index));
    }

    [[nodiscard]] inline auto
    node_indices_of(const GridCoordinates & cell_coordinates) const noexcept -> ElementNodes
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

    [[nodiscard]] inline auto
    node_indices_of(const CellIndex & index) const noexcept -> ElementNodes
    {
        return node_indices_of(cell_coordinates_at(index));
    }

    /** Get the number of distinct edges in this grid. **/
    [[nodiscard]] inline auto
    number_of_edges() const noexcept -> UInt
    {
        //        Grid         Slice           Row
        //      _ _ _ _       _ _ _ _
        //    /_/_/_/_/|     |_|_|_|_|      |_|_|_|_|
        //   /_/_/_/_/ |     |_|_|_|_|
        //   |_|_|_|_|/|     |_|_|_|_|
        //   |_|_|_|_|/|
        //   |_|_|_|_|/

        const auto & n = m_number_of_subdivisions;

        const auto & nx = n[0];
        const auto & ny = n[1];
        const auto & nz = n[2];

        const auto nb_edges_per_row = (2*nx + 1);
        const auto nb_edges_per_slice = (nb_edges_per_row * ny) + nx;
        const auto nb_edges_between_slices = (nx + 1)*(ny + 1);

        const auto nb_slices = (nz+1);

        const auto nb_edges_total = (nb_slices*nb_edges_per_slice) + (nz*nb_edges_between_slices);

        return nb_edges_total;
    }

    /** Get the node indices of the edge at edge index */
    [[nodiscard]] inline auto
    edge(const EdgeIndex & index) const noexcept -> EdgeNodes
    {
        //        Grid         Slice           Row
        //      _ _ _ _       _ _ _ _
        //    /_/_/_/_/|     |_|_|_|_|      |_|_|_|_|
        //   /_/_/_/_/ |     |_|_|_|_|
        //   |_|_|_|_|/|     |_|_|_|_|
        //   |_|_|_|_|/|
        //   |_|_|_|_|/

        const auto & n = m_number_of_subdivisions;

        const auto & nx = n[0];
        const auto & ny = n[1];

        const auto nb_edges_per_row = (2*nx + 1);
        const auto nb_edges_per_slice = (nb_edges_per_row * ny) + nx;
        const auto nb_edges_between_slices = (nx + 1)*(ny + 1);

        const auto slice = index / (nb_edges_per_slice + nb_edges_between_slices);
        const auto index_on_slice = index % (nb_edges_per_slice + nb_edges_between_slices);

        const bool edge_is_on_the_z_direction = index_on_slice > nb_edges_per_slice - 1;

        const auto slice_corner_node_index = slice * (nx+1) * (ny+1);

        NodeIndex n0, n1;

        if (edge_is_on_the_z_direction) {
            n0 = slice_corner_node_index + (index_on_slice - nb_edges_per_slice);
            n1 = n0 + nb_edges_between_slices;
        } else {
            const auto row = index_on_slice / nb_edges_per_row;
            const auto index_on_row = index_on_slice % nb_edges_per_row;
            const bool edge_is_on_the_x_direction = index_on_row < nx;

            if (edge_is_on_the_x_direction) {
                n0 = slice_corner_node_index + row * (nx+1) + index_on_row;
                n1 = n0 + 1;
            } else {
                n0 = slice_corner_node_index + row * (nx+1) + index_on_row - nx;
                n1 = n0 + (nx + 1);
            }
        }
        return {n0, n1};
    }

    /** Get the number of distinct faces in this grid. **/
    [[nodiscard]] inline auto
    number_of_faces() const noexcept -> UInt
    {
        //        Grid         Slice           Row
        //      _ _ _ _       _ _ _ _        _ _ _ _
        //    /_/_/_/_/|     |_|_|_|_|      |_|_|_|_|
        //   /_/_/_/_/ |     |_|_|_|_|
        //   |_|_|_|_|/|     |_|_|_|_|
        //   |_|_|_|_|/|
        //   |_|_|_|_|/

        const auto & n = Base::N();

        const auto & nx = n[0];
        const auto & ny = n[1];
        const auto & nz = n[2];

        const auto nb_faces_per_slice = nx*ny;
        const auto nb_faces_between_slices = ((nx + 1) * ny) + (nx * (ny + 1));

        const auto nb_slices = (nz+1);

        const auto nb_faces = (nb_slices*nb_faces_per_slice) + (nz*nb_faces_between_slices);

        return nb_faces;
    }

    /** Get the node indices of the face at the given face index */
    [[nodiscard]] inline auto
    face(const FaceIndex & index) const noexcept -> FaceNodes
    {
        //        Grid         Slice           Row
        //      _ _ _ _       _ _ _ _        _ _ _ _
        //    /_/_/_/_/|     |_|_|_|_|      |_|_|_|_|
        //   /_/_/_/_/ |     |_|_|_|_|
        //   |_|_|_|_|/|     |_|_|_|_|
        //   |_|_|_|_|/|
        //   |_|_|_|_|/

        const auto & n = Base::N();

        const auto & nx = n[0];
        const auto & ny = n[1];

        const auto nb_faces_per_row = nx;
        const auto nb_faces_per_slice = nx*ny;
        const auto nb_faces_between_slices = ((nx + 1) * ny) + (nx * (ny + 1));

        const auto slice = index / (nb_faces_per_slice + nb_faces_between_slices);
        const auto slice_corner_node_index = slice * (nx+1) * (ny+1);

        const auto index_on_slice = index % (nb_faces_per_slice + nb_faces_between_slices);
        const bool face_is_parallel_to_xy_plane = index_on_slice < nb_faces_per_slice;

        GridCoordinates e1; e1 << 1, 0, 0;
        GridCoordinates e2; e2 << 0, 1, 0;
        GridCoordinates e3; e3 << 0, 0, 1;

        if (face_is_parallel_to_xy_plane) {
            // (1) Faces on the slice (parallel to the x,y plane)
            // y
            // `
            // |4|5|6|7|
            // |0|1|2|3|...x

            const auto row = index_on_slice / nb_faces_per_row;
            const auto index_on_row = index_on_slice % nb_faces_per_row;

            NodeIndex n0 = slice_corner_node_index + row * (nx+1) + index_on_row;
            const auto corner_coordinates = node_coordinates_at(n0);
            return {{
                            node_index_at(corner_coordinates),
                            node_index_at(corner_coordinates + e1),
                            node_index_at(corner_coordinates + e1 + e2),
                            node_index_at(corner_coordinates + e2)
                    }};
        } else {
            const auto index_inbetween_slices = index_on_slice - nb_faces_per_slice; // 20 - 8 = 12
            const auto nb_of_faces_between_slices_xz = nb_faces_per_row * (ny+1); // 12
            const bool face_is_parallel_to_xz_plane = index_inbetween_slices < nb_of_faces_between_slices_xz;

            if (face_is_parallel_to_xz_plane) {
                // (2) Faces that are between two slices and that are parallel to the x,z plane
                const auto row = index_inbetween_slices / nb_faces_per_row;
                const auto index_on_row = index_inbetween_slices % nb_faces_per_row;
                const auto n0 = slice_corner_node_index + row * (nx+1) + index_on_row;
                const auto corner_coordinates = node_coordinates_at(n0);
                return {{
                    node_index_at(corner_coordinates),
                    node_index_at(corner_coordinates + e1),
                    node_index_at(corner_coordinates + e1 + e3),
                    node_index_at(corner_coordinates + e3)
                }};
            } else {
                // (3) Faces that are between two slices and that are parallel to the y,z plane
                const auto index_in_yz_faces = index_inbetween_slices - nb_of_faces_between_slices_xz;
                const auto n0 = slice_corner_node_index + index_in_yz_faces;
                const auto corner_coordinates = node_coordinates_at(n0);
                return {{
                    node_index_at(corner_coordinates),
                    node_index_at(corner_coordinates + e2),
                    node_index_at(corner_coordinates + e2 + e3),
                    node_index_at(corner_coordinates + e3)
                }};
            }
        }
    }
};

} // namespace caribou::topology

#endif //CARIBOU_TOPOLOGY_GRID_GRID_H
