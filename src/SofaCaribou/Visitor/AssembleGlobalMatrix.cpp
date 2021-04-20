#include <SofaCaribou/Visitor/AssembleGlobalMatrix.h>

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