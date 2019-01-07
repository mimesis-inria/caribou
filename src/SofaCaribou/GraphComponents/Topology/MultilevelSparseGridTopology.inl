#ifndef SOFACARIBOU_GRAPHCOMPONENTS_TOPOLOGY_MULTILEVELSPARSEGRIDTOPOLOGY_INL
#define SOFACARIBOU_GRAPHCOMPONENTS_TOPOLOGY_MULTILEVELSPARSEGRIDTOPOLOGY_INL

#include <SofaCaribou/GraphComponents/Topology/MultilevelSparseGridTopology.h>
#include <sofa/core/visual/VisualParams.h>

namespace SofaCaribou {
namespace GraphComponents {
namespace topology {

template <typename VecType>
MultilevelSparseGridTopology<VecType>::MultilevelSparseGridTopology()
        : DDGNode()
        , d_n(initData(&d_n, SofaVecInt(),"n","Grid resolution."))
        , d_number_of_subdivision(initData(&d_number_of_subdivision, "maximum_number_of_subdivision_levels", 0,
                                           "Number of subdivision levels of the boundary cells (one level split the cell in 8 subcells)."))
{
    addInput(&d_n);
    addInput(&d_number_of_subdivision);
}

template <typename VecType>
void MultilevelSparseGridTopology<VecType>::init() {
//    VecInt n = d_n.getValue();
//    Int maximum_number_of_subdivision_levels = d_number_of_subdivision.getValue();

    // Initialize the grid engine
//    p_grid.reset(
//        new GridType(
//            {0, 0, 0},       // anchor
//            n,               // subdivisions
//            {100, 100, 100}  // dimensions
//        )
//    );

//    const auto leaf_cells = p_grid->leaf_cells();
//    for (const auto & cell : leaf_cells) {
//
//    }
}

template <typename VecType>
void MultilevelSparseGridTopology<VecType>::onUpdate()
{
    // this doesn't work yet....
}


template <typename VecType>
void MultilevelSparseGridTopology<VecType>::draw(const sofa::core::visual::VisualParams* vparams)
{
    SOFA_UNUSED(vparams);
}

} // namespace topology
} // namespace GraphComponents
} // namespace SofaCaribou

#endif //SOFACARIBOU_GRAPHCOMPONENTS_TOPOLOGY_MULTILEVELSPARSEGRIDTOPOLOGY_INL
