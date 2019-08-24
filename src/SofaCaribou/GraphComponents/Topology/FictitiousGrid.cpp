#include <sofa/core/ObjectFactory.h>
#include "FictitiousGrid.inl"

#include <Caribou/Geometry/Triangle.h>

#include <iomanip>
#include <queue>

namespace SofaCaribou::GraphComponents::topology {

using sofa::defaulttype::Vec2Types;
using sofa::defaulttype::Vec3Types;

static const unsigned long long kelly_colors_hex[] = {
        0xFFFFB300, // Vivid Yellow
        0xFF803E75, // Strong Purple
        0xFFFF6800, // Vivid Orange
        0xFFA6BDD7, // Very Light Blue
        0xFFC10020, // Vivid Red
        0xFFCEA262, // Grayish Yellow
        0xFF817066, // Medium Gray

        // The following don't work well for people with defective color vision
        0xFF007D34, // Vivid Green
        0xFFF6768E, // Strong Purplish Pink
        0xFF00538A, // Strong Blue
        0xFFFF7A5C, // Strong Yellowish Pink
        0xFF53377A, // Strong Violet
        0xFFFF8E00, // Vivid Orange Yellow
        0xFFB32851, // Strong Purplish Red
        0xFFF4C800, // Vivid Greenish Yellow
        0xFF7F180D, // Strong Reddish Brown
        0xFF93AA00, // Vivid Yellowish Green
        0xFF593315, // Deep Yellowish Brown
        0xFFF13A13, // Vivid Reddish Orange
        0xFF232C16, // Dark Olive Green
};

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
        return;
    }

    // At this point, we have all the boundary cells, let's create clusters of cells regrouping neighbors cells of the same type
    begin = std::chrono::steady_clock::now();
    auto & regions = p_regions;
    regions.resize(0);
    auto & region_of_cell = p_region_of_cell;
    region_of_cell.resize(0);
    region_of_cell.resize(p_grid->number_of_cells(), -1);

    std::queue<CellIndex> remaining_cells;
    for (std::size_t i = 0; i < p_grid->number_of_cells(); ++i)
        remaining_cells.push(CellIndex(i));

    GridCoordinates upper_grid_boundary = p_grid->grid_coordinates_at(CellIndex(p_grid->number_of_cells()-1));

    const auto get_neighbors = [this, &upper_grid_boundary] (const GridCoordinates & coordinates) {
        std::array<CellIndex, 6> neighbors = {};
        for (std::size_t axis = 0; axis < 3; ++axis) {
            // Negative axis
            if (coordinates[axis] - 1 >= 0) {
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

        // Tag the current cell
        region_of_cell[cell_index] = current_region_index;
        regions[current_region_index].cells.emplace_back(cell_index);

        // Start the clustering algorithm
        std::queue<CellIndex> cluster_cells;
        cluster_cells.push(cell_index);

        while (not cluster_cells.empty()) {
            CellIndex cluster_cell_index = cluster_cells.front();
            cluster_cells.pop();

            // Add the neighbors that aren't already classified
            const auto neighbors = get_neighbors(p_grid->grid_coordinates_at(cluster_cell_index));
            for (const auto & neighbor_id : neighbors) {
                if (neighbor_id < 0)
                    continue;

                if (region_of_cell[neighbor_id] < 0 and regions[current_region_index].type == p_cells_types[neighbor_id]) {
                    // Tag the neighbor cell
                    region_of_cell[neighbor_id] = current_region_index;
                    regions[current_region_index].cells.emplace_back(neighbor_id);
                    cluster_cells.emplace(neighbor_id);
                }
            }
        }
    }

    end = std::chrono::steady_clock::now();
    msg_info() << "Computing the " << regions.size() << " cells regions in " << std::setprecision(3)
               << std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1000. / 1000.
               << " [ms]";

    // At this point, all cells have been classified into distinct regions of either boundary type or undefined type.

    // Next, the regions of undefined type which are surrounded by the grid's boundaries are tagged as outside cells
    begin = std::chrono::steady_clock::now();
    for (auto & region : p_regions) {
        if (region.type != Type::Undefined) {
            continue;
        }

        for (const auto & cell_id : region.cells) {
            const auto coordinates = p_grid->grid_coordinates_at(cell_id);
            bool cell_is_on_grid_boundary = false;
            for (std::size_t axis = 0; axis < 3; ++axis) {
                if (coordinates[axis] - 1 < 0 or coordinates[axis] + 1 > upper_grid_boundary[axis]) {
                    // We got a cell on the grid's boundary.
                    // We can tag all the region's cells as outside and exit the loop
                    cell_is_on_grid_boundary = true;
                    break;
                }
            }
            if (cell_is_on_grid_boundary) {
                for (const auto & internal_cell_id : region.cells) {
                    p_cells_types[internal_cell_id] = Type::Outside;
                }
                region.type = Type::Outside;
                break;
            }
        }
    }

    // Finally, we have both outside and boundary regions, the remaining regions are inside.
    for (auto & region : p_regions) {
        if (region.type != Type::Undefined) {
            continue;
        }

        region.type = Type::Inside;
        for (const auto & cell_id : region.cells) {
            p_cells_types[cell_id] = Type::Inside;
        }
    }
    end = std::chrono::steady_clock::now();
    msg_info() << "Computing the regions types in " << std::setprecision(3)
               << std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1000. / 1000.
               << " [ms]";

    static const std::string types_name[4] {
        "Undefined",
        "Inside",
        "Outside",
        "Boundary"
    };

    for (UNSIGNED_INTEGER_TYPE i = 0; i < regions.size(); ++i) {
        msg_info() << "Region #" << (i+1) << " has " << regions[i].cells.size()
                   << " cells of type " << types_name[(int)(regions[i].type)+1];
    }

}

template <>
void FictitiousGrid<Vec3Types>::draw(const sofa::core::visual::VisualParams* vparams)
{
    using Color = sofa::defaulttype::Vec4f;
    if (!p_grid)
        return;

    vparams->drawTool()->disableLighting();

    std::vector<sofa::defaulttype::Vector3> points_inside;
    std::vector<sofa::defaulttype::Vector3> points_boundary;
    for (UNSIGNED_INTEGER_TYPE i = 0; i < p_grid->number_of_nodes(); ++i) {
        const auto p = p_grid->node(i);
        if (p_node_types[i] == Type::Inside)
            points_inside.push_back({p[0], p[1], p[2]});
        else if (p_node_types[i] == Type::Boundary)
            points_boundary.push_back({p[0], p[1], p[2]});
    }
    vparams->drawTool()->drawPoints(points_inside, 5, Color(0, 1, 0, 1));
    vparams->drawTool()->drawPoints(points_boundary, 10, Color(1, 0, 0, 1));


    std::vector<sofa::defaulttype::Vector3> edge_points(p_grid->number_of_edges()*2);
    for (UNSIGNED_INTEGER_TYPE i = 0; i < p_grid->number_of_edges(); ++i) {
        const auto edge = p_grid->edge(i);
        edge_points[i*2] = {edge.node(0)[0], edge.node(0)[1], edge.node(0)[2]};
        edge_points[i*2+1] = {edge.node(1)[0], edge.node(1)[1], edge.node(1)[2]};
    }

    vparams->drawTool()->drawLines(edge_points, 0.5, Color(1,1,1, 1));

    if (! vparams->displayFlags().getShowWireFrame()) {
        for (UNSIGNED_INTEGER_TYPE i = 0; i < p_regions.size(); ++i) {
            const auto &region = p_regions[i];
            const auto &hex_color = kelly_colors_hex[i];
            const Color color(
                    (float) (((unsigned char)(hex_color >> 16)) / 255.),
                    (float) (((unsigned char)(hex_color >> 8)) / 255.),
                    (float) (((unsigned char)(hex_color >> 0)) / 255.),
                    (float) (0.75));

            std::vector<sofa::defaulttype::Vector3> nodes(region.cells.size() * 8);

            for (UNSIGNED_INTEGER_TYPE j = 0; j < region.cells.size(); ++j) {
                const auto &cell_index = region.cells[j];
                const auto cell = p_grid->cell_at(cell_index);
                for (UNSIGNED_INTEGER_TYPE k = 0; k < 8; ++k) {
                    nodes[j * 8 + k] = sofa::defaulttype::Vector3(cell.node(k)[0], cell.node(k)[1], cell.node(k)[2]);
                }
            }

            vparams->drawTool()->drawHexahedra(nodes, color);
        }
    }

    vparams->drawTool()->restoreLastState();
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
