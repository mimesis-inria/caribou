#pragma once

#include <sofa/simulation/MechanicalVisitor.h>

namespace SofaCaribou::visitor {

/**
 * Assemble the stiffness matrix (M+B+K) of every mechanical object into the mutli-matrix.
 *
 * Each top-level mechanical object (mo that aren't mapped to a higher level parent mo) has an entry
 * to its own matrix in the multi-matrix.
 *
 * This visitor will, in order:
 *   1. Assemble the stiffness matrix (M+B+K) of every top-level (not mapped) mechanical objects
 *      by calling `BaseForceField::addMBKToMatrix(mparams, matrix)` on every forcefields found
 *      in the context node of the mechanical object.
 *   2. Go down in the subtree of each top-level mechanical objects, and if a BaseMapping is
 *      found for which `BaseMapping::areMatricesMapped` is true, call
 *      `BaseForceField::addMBKToMatrix(mparams, matrix)` on every
 *      forcefields found in the context node of the mapped mechanical object.
 */
class AssembleGlobalMatrix : public sofa::simulation::MechanicalVisitor {
    using Base = sofa::simulation::MechanicalVisitor;
    using MechanicalParams = sofa::core::MechanicalParams;
    using MultiMatrixAccessor = sofa::core::behavior::MultiMatrixAccessor;
public:
    // Constructor
    AssembleGlobalMatrix(const MechanicalParams* mparams, const MultiMatrixAccessor* matrix )
    : Base(mparams), p_multi_matrix(matrix) {}

    Result fwdForceField(sofa::simulation::Node* node, sofa::core::behavior::BaseForceField* ff) override;

    bool stopAtMechanicalMapping(sofa::simulation::Node* node, sofa::core::BaseMapping* map) override;

    const char* getClassName() const override { return "AssembleGlobalMatrix"; }
private:
    const sofa::core::behavior::MultiMatrixAccessor * p_multi_matrix;
};

} // namespace SofaCaribou::visitor