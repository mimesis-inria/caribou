#include "ConstrainGlobalMatrix.h"

namespace SofaCaribou::visitor {
using namespace sofa::core;

auto ConstrainGlobalMatrix::fwdProjectiveConstraintSet(sofa::simulation::Node* /*node*/, sofa::core::behavior::BaseProjectiveConstraintSet* c) -> Result {
    c->applyConstraint(this->mparams, p_multi_matrix);
    return RESULT_CONTINUE;
}

bool ConstrainGlobalMatrix::stopAtMechanicalMapping(sofa::simulation::Node* /*node*/, sofa::core::BaseMapping* map) {
    return !map->areMatricesMapped();
}

} // namespace SofaCaribou::visitor