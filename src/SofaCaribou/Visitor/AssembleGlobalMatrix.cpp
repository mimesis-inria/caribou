#include <SofaCaribou/Visitor/AssembleGlobalMatrix.h>

DISABLE_ALL_WARNINGS_BEGIN
#include <sofa/core/behavior/BaseForceField.h>
DISABLE_ALL_WARNINGS_END

namespace SofaCaribou::visitor {
using namespace sofa::core;

auto AssembleGlobalMatrix::fwdForceField(sofa::simulation::Node* /*node*/, sofa::core::behavior::BaseForceField* ff) -> Result {
    ff->addMBKToMatrix(this->mparams, p_multi_matrix);
    return RESULT_CONTINUE;
}

bool AssembleGlobalMatrix::stopAtMechanicalMapping(sofa::simulation::Node* /*node*/, sofa::core::BaseMapping* map) {
    return !map->areMatricesMapped();
}

} // namespace SofaCaribou::visitor