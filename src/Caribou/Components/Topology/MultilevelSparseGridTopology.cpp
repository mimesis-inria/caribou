#include "MultilevelSparseGridTopology.h"
#include <sofa/core/ObjectFactory.h>

namespace caribou {
namespace components {
namespace topology {


} // namespace topology
} // namespace components
} // namespace caribou

SOFA_DECL_CLASS(MultilevelSparseGridTopology)
static int FEMForcefieldClass = sofa::core::RegisterObject("Caribou MultilevelSparseGridTopology")
                                        .add< caribou::components::topology::MultilevelSparseGridTopology >(true)
;