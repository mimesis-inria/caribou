#include <SofaCaribou/Ode/StaticODESolver.h>

#include <SofaCaribou/Visitor/AssembleGlobalMatrix.h>
#include <SofaCaribou/Visitor/ConstrainGlobalMatrix.h>

DISABLE_ALL_WARNINGS_BEGIN
#include <sofa/version.h>
#include <sofa/core/ConstraintParams.h>
#include <sofa/core/ObjectFactory.h>
#include <sofa/core/behavior/ConstraintSolver.h>
#include <sofa/helper/AdvancedTimer.h>
#include <sofa/simulation/MechanicalVisitor.h>
#if (defined(SOFA_VERSION) && SOFA_VERSION < 201299)
#include <sofa/simulation/MechanicalMatrixVisitor.h>
#else
#include <sofa/simulation/mechanicalvisitor/MechanicalApplyConstraintsVisitor.h>
#include <sofa/simulation/mechanicalvisitor/MechanicalComputeForceVisitor.h>
#include <sofa/simulation/mechanicalvisitor/MechanicalMultiVectorFromBaseVectorVisitor.h>
#include <sofa/simulation/mechanicalvisitor/MechanicalMultiVectorToBaseVectorVisitor.h>
#include <sofa/simulation/mechanicalvisitor/MechanicalPropagateOnlyPositionAndVelocityVisitor.h>
#include <sofa/simulation/mechanicalvisitor/MechanicalResetForceVisitor.h>
#include <sofa/simulation/mechanicalvisitor/MechanicalVOpVisitor.h>
using namespace sofa::simulation::mechanicalvisitor;
#endif
DISABLE_ALL_WARNINGS_END

namespace SofaCaribou::ode {

int StaticOdeSolverClass = sofa::core::RegisterObject("Static ODE Solver").add< StaticODESolver >();

using namespace sofa::simulation;
using Timer = sofa::helper::AdvancedTimer;

// Assemble F in A [dx] = F
void StaticODESolver::assemble_rhs_vector(const sofa::core::MechanicalParams &mechanical_parameters,
                                          const sofa::core::behavior::MultiMatrixAccessor & matrix_accessor,
                                          sofa::core::MultiVecDerivId & f_id,
                                          SofaCaribou::Algebra::BaseVector *f) {
    // 1. Clear the force vector (F := 0)
    MechanicalResetForceVisitor(&mechanical_parameters, f_id,
                                false /* onlyMapped */)
    .execute(this->getContext());

    // 2. Go down in the current context tree calling `addForce` on every force field components,
    //    then go up from the leaves calling `applyJT` on every mechanical mappings
    MechanicalComputeForceVisitor(&mechanical_parameters, f_id,
                                  true /* accumulate (to mapped node) */, true /*neglectingCompliance*/)
    .execute(this->getContext());

    // 3. Calls the "projectResponse" method of every `BaseProjectiveConstraintSet` objects found in the current
    //    context tree. For example, the `FixedConstraints` component will set entries of fixed nodes to zero.
    MechanicalApplyConstraintsVisitor(&mechanical_parameters, f_id,
                                      nullptr /*W (also project the given compliance matrix) */)
    .execute(this->getContext());

    // 4. Copy force vectors from every top level (unmapped) mechanical objects into the given system vector f
    MechanicalMultiVectorToBaseVectorVisitor(&mechanical_parameters, f_id /* source */, f /* destination */, &matrix_accessor)
    .execute(this->getContext());
}

// Assemble A in A [dx] = F
void StaticODESolver::assemble_system_matrix(const sofa::core::MechanicalParams &mechanical_parameters,
                                             sofa::core::behavior::DefaultMultiMatrixAccessor & matrix_accessor,
                                             SofaCaribou::Algebra::BaseMatrix *A) {
    // Step 1. Building stage:
    //         Here we go down on the current context sub-graph and call :
    //           1. ff->addMBKToMatrix(&K) for every force field "ff" found.
    //           2. pc->applyConstraint(&K) for every BaseProjectiveConstraintSet "pc" found.
    //         If a mechanical mapping "m" is found during the traversal, and m->areMatricesMapped() is false, the
    //         traversal stops in the subgraph of the mapping.
    matrix_accessor.setGlobalMatrix(A);
    sofa::core::MechanicalParams m_params (mechanical_parameters);
    m_params.setKFactor(-1.0);
    Timer::stepBegin("AssembleGlobalMatrix");
    visitor::AssembleGlobalMatrix(&m_params, &matrix_accessor).execute(this->getContext());
    Timer::stepEnd("AssembleGlobalMatrix");

    Timer::stepBegin("ConstrainGlobalMatrix");
    visitor::ConstrainGlobalMatrix(&m_params, &matrix_accessor).execute(this->getContext());
    Timer::stepEnd("ConstrainGlobalMatrix");

    // Step 2. Mechanical mappings
    //         In case we have mapped matrices, which is, system matrix of a slave mechanical object, accumulate its
    //         contribution to the global system matrix with:
    //           [A]ij += Jt * [A']ij * J
    //         where A is the master mechanical object's matrix, A' is the slave mechanical object matrix and J=m.getJ()
    //         is the mapping relation between the slave and its master.
    Timer::stepBegin("MappedMatrices");
    matrix_accessor.computeGlobalMatrix();
    Timer::stepEnd("MappedMatrices");

    // Step 3. Convert the system matrix to a compressed sparse matrix
    Timer::stepBegin("ConvertToSparse");
    A->compress();
    Timer::stepEnd("ConvertToSparse");

}

// Propagate Dx that was previously solved in A [dx] = F
void StaticODESolver::propagate_solution_increment(const sofa::core::MechanicalParams &mechanical_parameters,
                                                   const sofa::core::behavior::MultiMatrixAccessor & matrix_accessor,
                                                   const SofaCaribou::Algebra::BaseVector * dx,
                                                   sofa::core::MultiVecCoordId & x_id,
                                                   sofa::core::MultiVecDerivId & /* v_id */,
                                                   sofa::core::MultiVecDerivId & dx_id) {

    // 1. Copy vectors from the global system vector into every top level (unmapped) mechanical objects.
    MechanicalMultiVectorFromBaseVectorVisitor(&mechanical_parameters, dx_id, dx, &matrix_accessor).execute(this->getContext());

    // 2. x += dx
    MechanicalVOpVisitor(&mechanical_parameters, x_id, x_id, dx_id).execute(this->getContext());

    // 3. Calls "solveConstraint" method of every ConstraintSolver objects found in the current context tree.
    sofa::core::ConstraintParams constraint_parameters = mechanical_parameters;
    constraint_parameters.setOrder(sofa::core::ConstraintParams::POS);

    using Direction = sofa::core::objectmodel::BaseContext::SearchDirection;
    auto constraint_solvers = this->getContext()->getObjects<sofa::core::behavior::ConstraintSolver>(Direction::Local);
    for (auto * solver : constraint_solvers) {
        solver->solveConstraint(&constraint_parameters, x_id);
    }

    // 4. Propagate positions to mapped mechanical objects, for example, identity mappings, barycentric mappings, etc.
    //    This will call the methods apply and applyJ on every mechanical mappings.
    MechanicalPropagateOnlyPositionAndVelocityVisitor(&mechanical_parameters).execute(this->getContext());
}
} // namespace SofaCaribou::ode
