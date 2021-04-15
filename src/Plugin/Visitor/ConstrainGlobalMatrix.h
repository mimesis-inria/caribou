#pragma once

#include <SofaCaribou/config.h>
#include <sofa/simulation/MechanicalVisitor.h>

namespace SofaCaribou::visitor {

/**
 * Apply constraints on the stiffness matrix (M+B+K) of every mechanical object into the mutli-matrix.
 *
 * Each top-level mechanical object (mo that aren't mapped to a higher level parent mo) has an entry
 * to its own matrix in the multi-matrix.
 *
 * This visitor will, in order:
 *   1. Apply the constrains on the stiffness matrix (M+B+K) of every top-level (not mapped) mechanical
 *      objects by calling `BaseProjectiveConstraintSet::applyConstraint(mparams, matrix)` on every
 *      projective constraint sets found in the context node of the mechanical object.
 *   2. Go down in the subtree of each top-level mechanical objects, and if a BaseMapping is
 *      found for which `BaseMapping::areMatricesMapped` is true, call
 *      `BaseProjectiveConstraintSet::applyConstraint(mparams, matrix)` on every
 *      projective constraint sets found in the context node of the mapped mechanical object.
 */
class CARIBOU_API ConstrainGlobalMatrix : public sofa::simulation::MechanicalVisitor {
    using Base = sofa::simulation::MechanicalVisitor;
    using MechanicalParams = sofa::core::MechanicalParams;
    using MultiMatrixAccessor = sofa::core::behavior::MultiMatrixAccessor;
public:
    // Constructor
    ConstrainGlobalMatrix(const MechanicalParams* mparams, const MultiMatrixAccessor* matrix )
        : Base(mparams), p_multi_matrix(matrix) {}

    Result fwdProjectiveConstraintSet(sofa::simulation::Node * node, sofa::core::behavior::BaseProjectiveConstraintSet * c) override;

    bool stopAtMechanicalMapping(sofa::simulation::Node* node, sofa::core::BaseMapping* map) override;

    const char* getClassName() const override { return "ConstrainGlobalMatrix"; }
private:
    const sofa::core::behavior::MultiMatrixAccessor * p_multi_matrix;
};

} // namespace SofaCaribou::visitor