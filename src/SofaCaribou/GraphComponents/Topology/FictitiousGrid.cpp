#include <sofa/core/ObjectFactory.h>
#include "FictitiousGrid.inl"

#include <Caribou/Geometry/Triangle.h>

#include <iomanip>
#include <memory>
#include <queue>
#include <tuple>

namespace SofaCaribou::GraphComponents::topology {

//using sofa::defaulttype::Vec2Types;
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

//template<>
//void FictitiousGrid<Vec2Types>::create_grid()
//{
//    // Initializing the regular grid
//    const auto first_corner = d_min.getValue();
//    const auto second_corner = d_max.getValue();
//    const auto n = d_n.getValue();
//
//    WorldCoordinates anchor_position = {
//            std::min(first_corner[0], second_corner[0]),
//            std::min(first_corner[1], second_corner[1])
//    };
//    Dimensions grid_size = {
//            std::abs(first_corner[0] - second_corner[0]),
//            std::abs(first_corner[1] - second_corner[1])
//    };
//    Subdivisions grid_n = {n[0], n[1]};
//    p_grid = std::make_unique<GridType> (
//            anchor_position, grid_n, grid_size
//    );
//
//    p_node_types.resize(p_grid->number_of_nodes(), Type::Undefined);
//    p_cells_types.resize(p_grid->number_of_cells(), Type::Undefined);
//
//    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
//    if (d_use_implicit_surface.getValue() and p_implicit_test_callback) {
//        compute_cell_types_from_implicit_surface();
//    } else {
//        compute_cell_types_from_explicit_surface();
//    }
//    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
//    msg_info() << "Initialization of the grid finished in " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << " [ms]";
//}

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

    p_node_types.resize(p_grid->number_of_nodes(), Type::Undefined);
    p_cells_types.resize(p_grid->number_of_cells(), Type::Undefined);
    p_cells.resize(p_grid->number_of_cells());

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    if (d_use_implicit_surface.getValue() and p_implicit_test_callback) {
        compute_cell_types_from_implicit_surface();
    } else {
        compute_cell_types_from_explicit_surface();
    }
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    msg_info() << "Initialization of the grid finished in " << std::setprecision(3) << std::fixed
    << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << " [ms]";
}

//template<>
//void
//FictitiousGrid<Vec2Types>::compute_cell_types_from_explicit_surface()
//{
//    // We got a edge tesselation representation of the surface.
//    msg_error() << "Not yet implemented for 2D types.";
//}

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

    p_triangles_of_cell.resize(p_grid->number_of_cells());
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

        const caribou::geometry::Triangle<3> t(nodes[0], nodes[1], nodes[2]);

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
            const auto e = p_grid->cell_at(cell_index);
            begin = std::chrono::steady_clock::now();
            const bool intersects = e.intersects(t);
            end = std::chrono::steady_clock::now();
            time_to_find_intersections += std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count();
            if (intersects) {
                p_cells_types[cell_index] = Type::Boundary;
                p_triangles_of_cell[cell_index].emplace_back(triangle_index);
            }
        }
    }

    msg_info() << "Computing the bounding boxes of the surface elements in " << std::setprecision(3) << std::fixed
               << time_to_find_bounding_boxes/1000./1000. << " [ms]";
    msg_info() << "Computing the intersections with the surface in "  << std::setprecision(3) << std::fixed
               << time_to_find_intersections/1000./1000. << " [ms]";

    begin = std::chrono::steady_clock::now();
    subdivide_cells();
    end = std::chrono::steady_clock::now();
    msg_info() << "Computing the subdivisions in "  << std::setprecision(3) << std::fixed
               << std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count()/1000./1000. << " [ms]";

    if (!outside_triangles.empty()) {
        std::string triangle_indices = std::accumulate(std::next(outside_triangles.begin()), outside_triangles.end(),
            std::to_string(outside_triangles[0]),[](std::string s, const UNSIGNED_INTEGER_TYPE & index) {
            return std::move(s) + ", " + std::to_string(index);
        } );
        msg_error() << "Some triangles lie outside of the grid domain: " << triangle_indices;
        return;
    }

    // At this point, we have all the boundary cells, let's create regions of cells regrouping neighbors cells
    // of the same type. Special attention needs to be spent on boundary cells as if the number of subdivision is
    // greater than zero, we must also add subcells to the regions.
    begin = std::chrono::steady_clock::now();

    // First, add all leaf cells to a queue.
    std::queue<Cell*> remaining_cells;
    for (std::size_t i = 0; i < p_grid->number_of_cells(); ++i) {
        std::queue<Cell *> cells;
        cells.emplace(&p_cells[i]);
        while (not cells.empty()) {
            Cell * c = cells.front();
            cells.pop();

            if (c->childs) {
                for (Cell & child : (*c->childs)) {
                    cells.emplace(&child);
                }
            } else {
                remaining_cells.emplace(c);
            }

        }
    }

    // Iterate over every cells, and for each one, if it isn't yet classified,
    // start a clustering algorithm to fill the region's cells
    while (not remaining_cells.empty()) {
        Cell * c = remaining_cells.front();
        remaining_cells.pop();

        // If this cell was previously classified, skip it
        if (c->data->region_id > -1) {
            continue;
        }

        // Create a new region
        p_regions.push_back(Region {c->data->type, std::vector<Cell*> ()});
        std::size_t current_region_index = p_regions.size() - 1;

        // Tag the current cell
        c->data->region_id = current_region_index;
        p_regions[current_region_index].cells.push_back(c);

        // Start the clustering algorithm
        std::queue<Cell *> cluster_cells;
        cluster_cells.push(c);

        while (not cluster_cells.empty()) {
            Cell * cluster_cell = cluster_cells.front();
            cluster_cells.pop();

            // Add the neighbors that aren't already classified
            auto neighbors = get_neighbors(cluster_cell);
            for (Cell * neighbor_cell : neighbors) {
                if (neighbor_cell->data->region_id < 0 and
                    neighbor_cell->data->type == p_regions[current_region_index].type) {

                    // Tag the neighbor cell
                    neighbor_cell->data->region_id = current_region_index;
                    p_regions[current_region_index].cells.emplace_back(neighbor_cell);
                    cluster_cells.emplace(neighbor_cell);
                }
            }
        }
    }
    end = std::chrono::steady_clock::now();
    msg_info() << "Computing the " << p_regions.size() << " cells regions in " << std::fixed << std::setprecision(3)
               << std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1000. / 1000.
               << " [ms]";

    // At this point, all cells have been classified into distinct regions of either boundary type or undefined type.

    // Next, the regions of undefined type which are surrounded by the grid's boundaries are tagged as outside cells
    begin = std::chrono::steady_clock::now();
    for (auto & region : p_regions) {
        if (region.type != Type::Undefined) {
            continue;
        }

//        for (const auto & cell_id : region.cells) {
//            const auto coordinates = p_grid->grid_coordinates_at(cell_id);
//            bool cell_is_on_grid_boundary = false;
//            for (std::size_t axis = 0; axis < 3; ++axis) {
//                if (coordinates[axis] - 1 < 0 or coordinates[axis] + 1 > upper_grid_boundary[axis]) {
//                    // We got a cell on the grid's boundary.
//                    // We can tag all the region's cells as outside and exit the loop
//                    cell_is_on_grid_boundary = true;
//                    break;
//                }
//            }
//            if (cell_is_on_grid_boundary) {
//                for (const auto & internal_cell_id : region.cells) {
//                    p_cells_types[internal_cell_id] = Type::Outside;
//                }
//                region.type = Type::Outside;
//                break;
//            }
//        }
//    }
//
//    // Finally, we have both outside and boundary regions, the remaining regions are inside.
//    for (auto & region : p_regions) {
//        if (region.type != Type::Undefined) {
//            continue;
//        }
//
//        region.type = Type::Inside;
//        for (const auto & cell_id : region.cells) {
//            p_cells_types[cell_id] = Type::Inside;
//        }
    }
    end = std::chrono::steady_clock::now();
    msg_info() << "Computing the regions types in " << std::setprecision(3)
               << std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1000. / 1000.
               << " [ms]";
//
    begin = std::chrono::steady_clock::now();
    populate_drawing_vectors();
    end = std::chrono::steady_clock::now();
    msg_info() << "Populating the drawing vectors in " << std::setprecision(3) << std::fixed
               << std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1000. / 1000.
               << " [ms]";

//    static const std::string types_name[4] {
//        "undefined",
//        "inside",
//        "outside",
//        "boundary"
//    };
//
//    for (UNSIGNED_INTEGER_TYPE i = 0; i < regions.size(); ++i) {
//        msg_info() << "Region #" << (i+1) << " has " << regions[i].cells.size() << " "
//                   << types_name[(int)(regions[i].type)+1] << " cells.";
//    }
}

//template<>
//void
//FictitiousGrid<Vec2Types>::subdivide_cells()
//{
//    msg_error() << "Not yet implemented for 2D types.";
//}

template<>
void
FictitiousGrid<Vec3Types>::subdivide_cells()
{
    using Level = UNSIGNED_INTEGER_TYPE;
    using Weight = Float;

    const auto & number_of_subdivision = d_number_of_subdivision.getValue();
    const auto & surface_positions = d_surface_positions.getValue();
    const auto & surface_triangles = d_surface_triangles.getValue();

    for (UNSIGNED_INTEGER_TYPE cell_index = 0; cell_index < p_grid->number_of_cells(); ++cell_index) {
        const auto & triangles = p_triangles_of_cell[cell_index];
        std::queue<std::tuple<CellElement, Cell *, Weight, Level>> stack;

        // Initialize the stack with the current full cell
        Dimensions H;
        {
            p_cells[cell_index].index = cell_index;
            const CellElement cell = p_grid->cell_at(cell_index);
            const Weight weight = cell.jacobian().determinant();
            stack.emplace(p_grid->cell_at(cell_index), &p_cells[cell_index], weight, 0);
            H = cell.H();
        }

        while (not stack.empty()) {
            auto & s = stack.front();

            const CellElement & e = std::get<0>(s);
            Cell * cell = std::get<1>(s);
            const Weight & weight = std::get<2>(s);
            const Level & level = std::get<3>(s);

            Type type = Type::Undefined;
            bool subdivide_the_cell = false;

            // Checks if the current subcell intersects the boundary
            for (const auto & triangle_index : triangles) {
                const auto & triangle = surface_triangles[triangle_index];
                WorldCoordinates nodes [3];
                for (unsigned int i = 0; i < 3; ++i) {
                    const auto & node_index = triangle[i];

                    const Eigen::Map<const WorldCoordinates> p (&surface_positions[node_index][0]);
                    nodes[i] = p;
                }
                const caribou::geometry::Triangle<3> t(nodes[0], nodes[1], nodes[2]);
                if (e.intersects(t)) {
                    subdivide_the_cell = true;
                    type = Type::Boundary;
                    break;
                }
            }

            if (level+1 > number_of_subdivision or not subdivide_the_cell) {
                // We got a leaf, store the data
                cell->data = std::make_unique<CellData>(type, weight, -1);
            } else {
                // Split the cell into subcells
                const Weight w = weight / ((unsigned) 1<<Dimension);
                cell->childs = std::make_unique<std::array<Cell,(unsigned) 1 << Dimension>>();
                auto & childs = *(cell->childs);
                const auto & childs_elements = get_subcells(e);
                for (UNSIGNED_INTEGER_TYPE i = 0; i < childs.size(); ++i) {
                    childs[i].parent = cell;
                    childs[i].index = i;
                    stack.emplace(childs_elements[i], &(childs[i]), w, level+1);
                }
            }
            stack.pop();
        }
    }
}

//template<>
//std::array<FictitiousGrid<Vec2Types>::CellElement, (unsigned) 1 << FictitiousGrid<Vec2Types>::Dimension>
//FictitiousGrid<Vec2Types>::get_subcells(const CellElement & e) const
//{
//    msg_error() << "Not yet implemented for 2D types.";
//};

template<>
std::array<FictitiousGrid<Vec3Types>::CellElement, (unsigned) 1 << FictitiousGrid<Vec3Types>::Dimension>
FictitiousGrid<Vec3Types>::get_subcells(const CellElement & e) const
{
    const Dimensions h = e.H() /  2;
    return {{
        CellElement(e.T(LocalCoordinates {-0.5, -0.5, -0.5}), h),
        CellElement(e.T(LocalCoordinates {+0.5, -0.5, -0.5}), h),
        CellElement(e.T(LocalCoordinates {-0.5, +0.5, -0.5}), h),
        CellElement(e.T(LocalCoordinates {+0.5, +0.5, -0.5}), h),
        CellElement(e.T(LocalCoordinates {-0.5, -0.5, +0.5}), h),
        CellElement(e.T(LocalCoordinates {+0.5, -0.5, +0.5}), h),
        CellElement(e.T(LocalCoordinates {-0.5, +0.5, +0.5}), h),
        CellElement(e.T(LocalCoordinates {+0.5, +0.5, +0.5}), h)
    }};
};


template <>
void FictitiousGrid<Vec3Types>::draw(const sofa::core::visual::VisualParams* vparams)
{
    using Color = sofa::defaulttype::Vec4f;
    if (!p_grid)
        return;

    vparams->drawTool()->disableLighting();


//    vparams->drawTool()->drawPoints(p_drawing_nodes_vector, 1, Color(1, 0, 0, 1));
//    vparams->drawTool()->drawLines(p_drawing_edges_vector, 0.5, Color(1,1,1, 1));
    vparams->drawTool()->drawLines(p_drawing_subdivided_edges_vector, 1, Color(1,1,1, 1));

    if (! vparams->displayFlags().getShowWireFrame()) {
        for (UNSIGNED_INTEGER_TYPE i = 0; i < p_drawing_cells_vector.size(); ++i) {
            const auto &hex_color = kelly_colors_hex[i];
            const Color color(
                    (float) (((unsigned char)(hex_color >> (unsigned) 16)) / 255.),
                    (float) (((unsigned char)(hex_color >> (unsigned) 8)) / 255.),
                    (float) (((unsigned char)(hex_color >> (unsigned) 0)) / 255.),
                    (float) (0.75));

            vparams->drawTool()->drawHexahedra(p_drawing_cells_vector[i], color);
        }
    }

    vparams->drawTool()->restoreLastState();
}

// This will force the compiler to compile the class with some template type
//template class FictitiousGrid<Vec2Types>;
template class FictitiousGrid<Vec3Types>;

// Add the sofa component to the object factory
int FictitiousGridClass = sofa::core::RegisterObject("Caribou FictitiousGrid")
//        .add< FictitiousGrid<Vec2Types> >()
        .add< FictitiousGrid<Vec3Types> >(true)
;

} // namespace SofaCaribou::GraphComponents::topology
