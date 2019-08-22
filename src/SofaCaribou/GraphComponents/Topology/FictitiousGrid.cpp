#include <sofa/core/ObjectFactory.h>
#include "FictitiousGrid.inl"

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
    for (const auto & triangle : triangles) {
        WorldCoordinates nodes [3];
        for (unsigned int i = 0; i < 3; ++i) {
            const auto & node_index = triangle[i];
            if (node_index >= positions.size()) {
                msg_error() << "Some triangles have their node index greater than the size of the position vector.";
                return;
            }

            const Eigen::Map<const WorldCoordinates> p (&positions[node_index][0]);
            const auto cell_index = p_grid->cell_index_at(WorldCoordinates(p));
            if (cell_index < 0 || cell_index > (signed) p_grid->number_of_cells()) {
                msg_error() << "Some triangles lie outside of the grid domain.";
                return;
            }

            nodes[i] = p;
        }

        // Get all the cells enclosing the three nodes of the triangles
        const auto enclosing_cells = p_grid->cells_enclosing(nodes[0], nodes[1], nodes[2]);
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
