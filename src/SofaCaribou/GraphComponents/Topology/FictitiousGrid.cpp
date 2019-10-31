#include <sofa/core/ObjectFactory.h>

#ifdef CARIBOU_WITH_OPENMP
#include <omp.h>
#endif

#include "FictitiousGrid.inl"

#include <Caribou/Geometry/Triangle.h>

namespace SofaCaribou::GraphComponents::topology {

//using sofa::defaulttype::Vec2Types;
using sofa::defaulttype::Vec3Types;

template<>
const FictitiousGrid<Vec2Types>::GridCoordinates FictitiousGrid<Vec2Types>::subcell_coordinates[]  = {
        {0, 0},
        {1, 0},
        {0, 1},
        {1, 1}
};

template<>
const FictitiousGrid<Vec3Types>::GridCoordinates FictitiousGrid<Vec3Types>::subcell_coordinates[]  = {
        {0, 0, 0},
        {1, 0, 0},
        {0, 1, 0},
        {1, 1, 0},
        {0, 0, 1},
        {1, 0, 1},
        {0, 1, 1},
        {1, 1, 1},
};

template<>
void
FictitiousGrid<Vec2Types>::tag_intersected_cells()
{
    msg_error() << "Not yet implemented for 2D.";
}

template<>
void
FictitiousGrid<Vec3Types>::tag_intersected_cells()
{
    BEGIN_CLOCK;
    // We got a triangle tesselation representation of the surface.
    const auto & positions = d_surface_positions.getValue();
    const auto & triangles = d_surface_triangles.getValue();
    int64_t time_to_find_bounding_boxes = 0;
    int64_t time_to_find_intersections = 0;

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
        TICK;
        const auto enclosing_cells = p_grid->cells_enclosing(nodes[0], nodes[1], nodes[2]);
        time_to_find_bounding_boxes += TOCK;

        if (enclosing_cells.empty()) {
            msg_error() << "Triangle #"<< triangle_index << " has no enclosing cells.";
            return;
        }
        if (enclosing_cells.size() == 1) {
            p_cells_types[*enclosing_cells.begin()] = Type::Boundary;
            p_triangles_of_cell[*enclosing_cells.begin()].emplace_back(triangle_index);
        } else {
            for (const auto &cell_index : enclosing_cells) {
                const auto e = p_grid->cell_at(cell_index);
                TICK;
                const bool intersects = e.intersects(t);
                time_to_find_intersections += TOCK;
                if (intersects) {
                    p_cells_types[cell_index] = Type::Boundary;
                    p_triangles_of_cell[cell_index].emplace_back(triangle_index);
                }
            }
        }
    }

    msg_info() << "Computing the bounding boxes of the surface elements in " << std::setprecision(3) << std::fixed
               << time_to_find_bounding_boxes/1000./1000. << " [ms]";
    msg_info() << "Computing the intersections with the surface in "  << std::setprecision(3) << std::fixed
               << time_to_find_intersections/1000./1000. << " [ms]";

    if (!outside_triangles.empty()) {
        std::string triangle_indices = std::accumulate(std::next(outside_triangles.begin()), outside_triangles.end(),
                                                       std::to_string(outside_triangles[0]),[](std::string s, const UNSIGNED_INTEGER_TYPE & index) {
                return std::move(s) + ", " + std::to_string(index);
            } );
        msg_error() << "Some triangles lie outside of the grid domain: " << triangle_indices;
        return;
    }
}

template<>
void
FictitiousGrid<Vec2Types>::subdivide_intersected_cells()
{
    BEGIN_CLOCK;
    using Level = UNSIGNED_INTEGER_TYPE;
    using Weight = Float;

    const auto & number_of_subdivision = d_number_of_subdivision.getValue();
    bool use_implicit_surface = (d_use_implicit_surface.getValue() and p_implicit_test_callback);

    if (not use_implicit_surface) {
        msg_error() << "Not yet implemented for 2D types.";
        return;
    }

    TICK;
    for (UNSIGNED_INTEGER_TYPE cell_index = 0; cell_index < p_grid->number_of_cells(); ++cell_index) {
        std::queue<std::tuple<CellElement, Cell *, Weight, Level>> stack;

        // Initialize the stack with the current full cell
        Dimensions H;
        {
            p_cells[cell_index].index = cell_index;
            const CellElement cell = p_grid->cell_at(cell_index);
            const Weight weight = 1;
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

            if (use_implicit_surface) {
                constexpr UNSIGNED_INTEGER_TYPE INSIDE = 0;
                constexpr UNSIGNED_INTEGER_TYPE OUTSIDE = 1;
                constexpr UNSIGNED_INTEGER_TYPE BOUNDARY = 2;
                UNSIGNED_INTEGER_TYPE types[3] = {0, 0, 0};
                for (UNSIGNED_INTEGER_TYPE i = 0; i < caribou::traits<CellElement>::NumberOfNodes; ++i) {
                    const auto t = p_implicit_test_callback(e.node(i));
                    if (t < 0)
                        types[INSIDE]++;
                    else if (t > 0)
                        types[OUTSIDE]++;
                    else
                        types[BOUNDARY]++;
                }

                if (types[INSIDE] == caribou::traits<CellElement>::NumberOfNodes) {
                    type = Type::Inside;
                } else if (types[OUTSIDE] == caribou::traits<CellElement>::NumberOfNodes) {
                    type = Type::Outside;
                } else {
                    type = Type::Boundary;
                    subdivide_the_cell = true;
                }
            } else {
                // Todo: intersections quads-edges
            }

            if (level+1 > number_of_subdivision or not subdivide_the_cell) {
                // We got a leaf, store the data
                cell->data = std::make_unique<CellData>(type, weight, -1);
            } else {
                // Split the cell into subcells
                cell->data.reset();
                const Weight w = weight / ((unsigned) 1<<Dimension);
                cell->childs = std::make_unique<std::array<Cell,(unsigned) 1 << Dimension>>();
                auto & childs = *(cell->childs);
                const auto & childs_elements = get_subcells_elements(e);
                for (UNSIGNED_INTEGER_TYPE i = 0; i < childs.size(); ++i) {
                    childs[i].parent = cell;
                    childs[i].index = i;
                    stack.emplace(childs_elements[i], &(childs[i]), w, level+1);
                }
            }
            stack.pop();
        }
    }
    msg_info() << "Computing the subdivisions in "  << std::setprecision(3) << std::fixed
               << TOCK/1000./1000. << " [ms]";
}

template<>
void
FictitiousGrid<Vec3Types>::subdivide_intersected_cells()
{
    BEGIN_CLOCK;
    using Level = UNSIGNED_INTEGER_TYPE;
    using Weight = Float;

    const auto & number_of_subdivision = d_number_of_subdivision.getValue();
    const auto & surface_positions = d_surface_positions.getValue();
    const auto & surface_triangles = d_surface_triangles.getValue();
    bool use_implicit_surface = (d_use_implicit_surface.getValue() and p_implicit_test_callback);

    TICK;
#pragma omp parallel for default(none) shared(use_implicit_surface, surface_triangles, surface_positions, number_of_subdivision)
    for (UNSIGNED_INTEGER_TYPE cell_index = 0; cell_index < p_grid->number_of_cells(); ++cell_index) {
        const auto & triangles = p_triangles_of_cell[cell_index];
        std::queue<std::tuple<CellElement, Cell *, Weight, Level>> stack;

        // Initialize the stack with the current full cell
        Dimensions H;
        {
            p_cells[cell_index].index = cell_index;
            const CellElement cell = p_grid->cell_at(cell_index);
            const Weight weight = 1;
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

            if (use_implicit_surface) {
                constexpr UNSIGNED_INTEGER_TYPE INSIDE = 0;
                constexpr UNSIGNED_INTEGER_TYPE OUTSIDE = 1;
                constexpr UNSIGNED_INTEGER_TYPE BOUNDARY = 2;

                UNSIGNED_INTEGER_TYPE types[3] = {0, 0, 0};
                for (UNSIGNED_INTEGER_TYPE i = 0; i < caribou::traits<CellElement>::NumberOfNodes; ++i) {
                    const auto t = p_implicit_test_callback(e.node(i));
                    if (t < 0)
                        types[INSIDE]++;
                    else if (t > 0)
                        types[OUTSIDE]++;
                    else
                        types[BOUNDARY]++;
                }

                if (types[INSIDE] == caribou::traits<CellElement>::NumberOfNodes) {
                    type = Type::Inside;
                } else if (types[OUTSIDE] == caribou::traits<CellElement>::NumberOfNodes) {
                    type = Type::Outside;
                } else {
                    type = Type::Boundary;
                    subdivide_the_cell = true;
                }
            } else {
                for (const auto &triangle_index : triangles) {
                    const auto &triangle = surface_triangles[triangle_index];
                    WorldCoordinates nodes[3];
                    for (unsigned int i = 0; i < 3; ++i) {
                        const auto &node_index = triangle[i];

                        const Eigen::Map<const WorldCoordinates> p(&surface_positions[node_index][0]);
                        nodes[i] = p;
                    }
                    const caribou::geometry::Triangle<3> t(nodes[0], nodes[1], nodes[2]);
                    const bool intersects = e.intersects(t);

                    if (intersects) {
                        subdivide_the_cell = true;
                        type = Type::Boundary;
                        break;
                    }
                }
            }

            if (level+1 > number_of_subdivision or not subdivide_the_cell) {
                // We got a leaf, store the data
                cell->data = std::make_unique<CellData>(type, weight, -1);
            } else {
                // Split the cell into subcells
                cell->data.reset();
                const Weight w = weight / ((unsigned) 1<<Dimension);
                cell->childs = std::make_unique<std::array<Cell,(unsigned) 1 << Dimension>>();
                auto & childs = *(cell->childs);
                const auto & childs_elements = get_subcells_elements(e);
                for (UNSIGNED_INTEGER_TYPE i = 0; i < childs.size(); ++i) {
                    childs[i].parent = cell;
                    childs[i].index = i;
                    stack.emplace(childs_elements[i], &(childs[i]), w, level+1);
                }
            }
            stack.pop();
        }
    }
    msg_info() << "Computing the subdivisions in "  << std::setprecision(3) << std::fixed
               << TOCK/1000./1000. << " [ms]";
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
