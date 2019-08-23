#include <sofa/core/ObjectFactory.h>
#include "FictitiousGrid.inl"

#include <Caribou/Geometry/Triangle.h>

#include <iomanip>
#include <queue>

namespace SofaCaribou::GraphComponents::topology {

using sofa::defaulttype::Vec2Types;
using sofa::defaulttype::Vec3Types;

template<>
void FictitiousGrid<Vec2Types>::create_grid()
{
    // Initializing the regular grid
    const auto first_corner = d_min.getValue();
    const auto second_corner = d_max.getValue();
    const auto n = d_n.getValue();

    WorldCoordinates anchor_position = {
            std::min(first_corner[0], second_corner[0]),
            std::min(first_corner[1], second_corner[1])
    };
    Dimensions grid_size = {
            std::abs(first_corner[0] - second_corner[0]),
            std::abs(first_corner[1] - second_corner[1])
    };
    Subdivisions grid_n = {n[0], n[1]};
    p_grid = std::make_unique<GridType> (
            anchor_position, grid_n, grid_size
    );
}

template<>
void FictitiousGrid<Vec3Types>::create_grid()
{
    // Initializing the regular grid
    const auto first_corner = d_min.getValue();
    const auto second_corner = d_max.getValue();
    const auto n = d_n.getValue();

    WorldCoordinates anchor_position = {
            std::min(first_corner[0], second_corner[0]),
            std::min(first_corner[1], second_corner[1]),
            std::min(first_corner[2], second_corner[2])
    };
    Dimensions grid_size = {
            std::abs(first_corner[0] - second_corner[0]),
            std::abs(first_corner[1] - second_corner[1]),
            std::abs(first_corner[2] - second_corner[2])
    };
    Subdivisions grid_n = {n[0], n[1], n[2]};
    p_grid = std::make_unique<GridType> (
            anchor_position, grid_n, grid_size
    );
}

template<>
void
FictitiousGrid<Vec2Types>::compute_cell_types_from_explicit_surface()
{
    // We got a edge tesselation representation of the surface.
    msg_error() << "Not yet implemented for 2D types.";
}

template<>
void
FictitiousGrid<Vec3Types>::compute_cell_types_from_explicit_surface()
{
    // We got a triangle tesselation representation of the surface.
    const auto & positions = d_surface_positions.getValue();
    const auto & triangles = d_surface_triangles.getValue();
    int64_t time_to_find_bounding_boxes = 0;
    int64_t time_to_find_intersections = 0;
    std::chrono::steady_clock::time_point begin, end;

    std::vector<UNSIGNED_INTEGER_TYPE> outside_triangles;
    for (std::size_t triangle_index = 0; triangle_index < triangles.size(); ++triangle_index) {
        const auto & triangle = triangles[triangle_index];
        WorldCoordinates nodes [3];
        bool triangle_is_outside = false;
        for (unsigned int i = 0; i < 3; ++i) {
            const auto & node_index = triangle[i];
            if (node_index >= positions.size()) {
                msg_error() << "Some triangles have their node index greater than the size of the position vector.";
                return;
            }

            const Eigen::Map<const WorldCoordinates> p (&positions[node_index][0]);
            if (!p_grid->contains(p)) {
                triangle_is_outside = true;
            }

            nodes[i] = p;
        }

        if (triangle_is_outside) {
            outside_triangles.push_back(triangle_index);
            continue;
        }

        caribou::geometry::Triangle<3> t(nodes[0], nodes[1], nodes[2]);

        // Get all the cells enclosing the three nodes of the triangles
        begin = std::chrono::steady_clock::now();
        const auto enclosing_cells = p_grid->cells_enclosing(nodes[0], nodes[1], nodes[2]);
        end = std::chrono::steady_clock::now();
        time_to_find_bounding_boxes += std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count();

        if (enclosing_cells.empty()) {
            msg_error() << "Triangle #"<< triangle_index << " has no enclosing cells.";
            return;
        }
        for (const auto & cell_index : enclosing_cells) {
            const auto cell = p_grid->cell_at(cell_index);
            begin = std::chrono::steady_clock::now();
            const bool intersects = cell.intersects(t);
            end = std::chrono::steady_clock::now();
            time_to_find_intersections += std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count();
            if (intersects) {
                p_cells_types[cell_index] = Type::Boundary;

                const auto nodes_indices = p_grid->node_indices_of(cell_index);
                for (const auto & node_index : nodes_indices) {
                    p_node_types[node_index] = Type::Boundary;
                }
            }
        }
    }

    msg_info() << "Computing the bounding boxes of the surface elements in " << std::setprecision(3) << time_to_find_bounding_boxes/1000/1000 << " [ms]";
    msg_info() << "Computing the intersections with the surface in "  << std::setprecision(3) << time_to_find_intersections/1000/1000 << " [ms]";

    if (!outside_triangles.empty()) {
        std::string triangle_indices = std::accumulate(std::next(outside_triangles.begin()), outside_triangles.end(), std::to_string(outside_triangles[0]),[](std::string s, const UNSIGNED_INTEGER_TYPE & index) {
            return std::move(s) + ", " + std::to_string(index);
        } );
        msg_error() << "Some triangles lie outside of the grid domain: " << triangle_indices;
    }

    // At this point, we have all the boundary cells, let's create clusters of cells regrouping neighbors cells of the same type
    struct Region {
        Type type = Type::Undefined;
        std::vector<CellIndex> cells;
    };

    std::vector<Region> regions;
    std::vector<int> region_of_cell (p_grid->number_of_cells(), -1);

    std::queue<CellIndex> remaining_cells;
    for (std::size_t i = 0; i < p_grid->number_of_cells(); ++i)
        remaining_cells.push(CellIndex(i));

    GridCoordinates upper_grid_boundary = p_grid->grid_coordinates_at(CellIndex(p_grid->number_of_cells()-1));

    const auto get_neighbors = [this, &upper_grid_boundary] (const GridCoordinates & coordinates) {
        std::array<CellIndex, 6> neighbors = {};
        for (std::size_t axis = 0; axis < 3; ++axis) {
            // Negative axis
            if (coordinates[axis] - 1 > 0) {
                GridCoordinates c = coordinates;
                c[axis] -= 1;
                neighbors[axis * 2] = p_grid->cell_index_at(c);
            } else {
                neighbors[axis * 2] = -1;
            }

            // Positive axis
            if (coordinates[axis] + 1 <= upper_grid_boundary[axis]) {
                GridCoordinates c = coordinates;
                c[axis] += 1;
                neighbors[axis * 2 + 1] = p_grid->cell_index_at(c);
            } else {
                neighbors[axis * 2 + 1] = -1;
            }
        }

        return neighbors;
    };

    // Iterate over every cells, and for each one, if it isn't yet classified,
    // start a clustering algorithm to fill the region's cells
    while (not remaining_cells.empty()) {
        CellIndex cell_index = remaining_cells.front();
        remaining_cells.pop();

        // If this cell was previously classified, skip it
        if (region_of_cell[cell_index] > -1) {
            continue;
        }

        // Create a new region
        regions.push_back(Region {p_cells_types[cell_index], std::vector<CellIndex> ()});
        std::size_t current_region_index = regions.size() - 1;

        // Start the clustering algorithm
        std::queue<CellIndex> cluster_cells;
        cluster_cells.push(cell_index);

        while (not cluster_cells.empty()) {
            CellIndex cluster_cell_index = cluster_cells.front();
            cluster_cells.pop();

            // Tag the current cell
            region_of_cell[cluster_cell_index] = current_region_index;
            regions[current_region_index].cells.emplace_back(cluster_cell_index);

            // Add the neighbors that aren't already classified
            const auto neighbors = get_neighbors(p_grid->grid_coordinates_at(cluster_cell_index));
            for (const auto & neighbor_id : neighbors) {
                if (neighbor_id < 0)
                    continue;

                if (region_of_cell[neighbor_id] < 0 and regions[current_region_index].type == p_cells_types[neighbor_id]) {
                    cluster_cells.emplace(neighbor_id);
                }
            }
        }
    }

    // At this point, all cells have been classified into distinct regions
    msg_info() << "The grid contains " << regions.size() << " regions:";
    for (std::size_t i = 0; i < regions.size(); ++i) {
        msg_info() << "Region #" << i+1 << " has " << regions[i].cells.size() << " cells";
    }
}

// This will force the compiler to compile the class with some template type
template class FictitiousGrid<Vec2Types>;
template class FictitiousGrid<Vec3Types>;

// Add the sofa component to the object factory
int FictitiousGridClass = sofa::core::RegisterObject("Caribou FictitiousGrid")
        .add< FictitiousGrid<Vec2Types> >()
        .add< FictitiousGrid<Vec3Types> >(true)
;

} // namespace SofaCaribou::GraphComponents::topology
