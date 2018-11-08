#include "MultilevelSparseGridTopology.h"
#include <sofa/core/ObjectFactory.h>

namespace caribou {
namespace components {
namespace topology {


MultilevelSparseGridTopology::MultilevelSparseGridTopology()
        : Parent(),
        d_number_of_subdivision(initData(&d_number_of_subdivision, "number_of_subdivision_levels", 0,
                                         "Number of subdivision levels of the boundary cells (one level split the cell in 8 subcells).")) {
}

void MultilevelSparseGridTopology::init() {
    Parent::init();

    m_cells.resize(getNbHexahedra());

    for (HexaID hexa_id = 0; hexa_id < getNbHexahedra(); ++hexa_id) {
    }


}

} // namespace topology
} // namespace components
} // namespace caribou

SOFA_DECL_CLASS(MultilevelSparseGridTopology)
static int FEMForcefieldClass = sofa::core::RegisterObject("Caribou MultilevelSparseGridTopology")
                                        .add< caribou::components::topology::MultilevelSparseGridTopology >(true)
;