#ifndef SOFACARIBOU_GRAPHCOMPONENTS_TOPOLOGY_FICTITIOUSGRID_INL
#define SOFACARIBOU_GRAPHCOMPONENTS_TOPOLOGY_FICTITIOUSGRID_INL

#include <SofaCaribou/GraphComponents/Topology/FictitiousGrid.h>
#include <sofa/core/visual/VisualParams.h>

namespace SofaCaribou::GraphComponents::topology {

template <typename VecType>
FictitiousGrid<VecType>::FictitiousGrid()
        : d_n(initData(&d_n, SofaVecInt(),"n","Grid resolution."))
        , d_min(initData(&d_min, SofaVecFloat(),"min","First corner node position of the grid's bounding box."))
        , d_max(initData(&d_max, SofaVecFloat(),"max","Second corner node position of the grid's bounding box."))
        , d_number_of_subdivision(initData(&d_number_of_subdivision, "maximum_number_of_subdivision_levels", 0,
                                           "Number of subdivision levels of the boundary cells (one level split the cell in 8 subcells)."))
{
}

template <typename VecType>
void FictitiousGrid<VecType>::init() {

    // Initializing the regular grid
    const auto first_corner = d_min.getValue();
    const auto second_corner = d_min.getValue();
    const auto n = d_n.getValue();

    if constexpr (Dimension == 2) {
        WorldCoordinates anchor_position = {
                std::min(first_corner[0], second_corner[0]),
                std::min(first_corner[1], second_corner[1])
        };
        Dimensions grid_size = {
                std::abs(first_corner[0] - second_corner[0]),
                std::abs(first_corner[1] - second_corner[1])
        };
        Subdivisions grid_n = {n[0], n[1]};
        p_grid.reset(
                new GridType(anchor_position, grid_n, grid_size)
        );
    } else { // Dimension == 3
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
        p_grid.reset(
                new GridType(anchor_position, grid_n, grid_size)
        );

    }

    // Set all cells at an undefined state (will be defined later through seeding and intersections with a boundary mesh)
    p_cells.resize(p_grid->number_of_cells());
}

template <typename VecType>
void FictitiousGrid<VecType>::draw(const sofa::core::visual::VisualParams* vparams)
{
    SOFA_UNUSED(vparams);
}

} // namespace SofaCaribou::GraphComponents::topology

#endif //SOFACARIBOU_GRAPHCOMPONENTS_TOPOLOGY_FICTITIOUSGRID_INL
