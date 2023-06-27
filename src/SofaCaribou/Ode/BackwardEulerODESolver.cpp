#include <SofaCaribou/Ode/BackwardEulerODESolver.h>

#include <SofaCaribou/Visitor/AssembleGlobalMatrix.h>
#include <SofaCaribou/Visitor/ConstrainGlobalMatrix.h>

DISABLE_ALL_WARNINGS_BEGIN
#include <sofa/version.h>
#include <sofa/core/ConstraintParams.h>
#include <sofa/core/ObjectFactory.h>
#include <sofa/core/behavior/ConstraintSolver.h>
#include <sofa/helper/AdvancedTimer.h>
#include <sofa/simulation/MechanicalVisitor.h>
#include <sofa/simulation/VectorOperations.h>
#if (defined(SOFA_VERSION) && SOFA_VERSION < 201299)
#include <sofa/simulation/MechanicalMatrixVisitor.h>
#else
#include <sofa/simulation/mechanicalvisitor/MechanicalAddMBKdxVisitor.h>
#include <sofa/simulation/mechanicalvisitor/MechanicalApplyConstraintsVisitor.h>
#include <sofa/simulation/mechanicalvisitor/MechanicalComputeForceVisitor.h>
#include <sofa/simulation/mechanicalvisitor/MechanicalMultiVectorFromBaseVectorVisitor.h>
#include <sofa/simulation/mechanicalvisitor/MechanicalMultiVectorToBaseVectorVisitor.h>
#include <sofa/simulation/mechanicalvisitor/MechanicalPropagateOnlyPositionAndVelocityVisitor.h>
#include <sofa/simulation/mechanicalvisitor/MechanicalResetForceVisitor.h>
#include <sofa/simulation/mechanicalvisitor/MechanicalVOpVisitor.h>
using namespace sofa::simulation::mechanicalvisitor;
#endif
DISABLE_ALL_WARNINGS_BEGIN

namespace SofaCaribou::ode {

int BackwardEulerClass = sofa::core::RegisterObject("Backward Euler ODE Solver").add< BackwardEulerODESolver >();

using namespace sofa::simulation;
using sofa::core::behavior::MultiMatrixAccessor;
using sofa::core::MechanicalParams;
using sofa::core::MultiVecCoordId;
using sofa::core::MultiVecDerivId;
using sofa::core::behavior::DefaultMultiMatrixAccessor;
using Timer = sofa::helper::AdvancedTimer;

// Constructor
BackwardEulerODESolver::BackwardEulerODESolver()
: d_rayleigh_stiffness(initData(&d_rayleigh_stiffness,
    (double) 0.0,
    "rayleigh_stiffness",
    "The stiffness factor 'r_k' used in the Rayleigh's damping matrix `D = r_m M + r_k K`."))
, d_rayleigh_mass(initData(&d_rayleigh_mass,
    (double) 0.0,
    "rayleigh_mass",
    "The mass factor 'r_m' used in the Rayleigh's damping matrix `D = r_m M + r_k K`."))
{}


void BackwardEulerODESolver::solve(const sofa::core::ExecParams *params, SReal dt, sofa::core::MultiVecCoordId x_id,
                                   sofa::core::MultiVecDerivId v_id) {
    // Save up the current position and velocity multi vectors in order to reuse it during the time stepping part
    sofa::core::MechanicalParams mechanical_parameters (*params);
    sofa::simulation::common::VectorOperations vop( &mechanical_parameters, this->getContext() );
    vop.v_realloc(p_previous_x_id, false /* interactionForceField */, true /* propagate [to mapped MO] */);
    vop.v_eq(p_previous_x_id, x_id); // x_0 = x

    vop.v_realloc(p_previous_v_id, false /* interactionForceField */, true /* propagate [to mapped MO] */);
    vop.v_eq(p_previous_v_id, v_id); // v_0 = v

    // Allocate the acceleration vector
    vop.v_realloc(p_a_id, false /* interactionForceField */, true /* propagate [to mapped MO] */);
    vop.v_clear(p_a_id);

    // Let the NR do its job
    NewtonRaphsonSolver::solve(params, dt, x_id, v_id);
}


// Assemble F in A [da] = F
//          F =   [(1 + h r_m) M  +  h C  +  h r_k K ] a
//              + [r_m M  +  C  +  r_k K] v
//              + [R(x_{n+1}) - P_n]
//            = f2 + f1 + f0
void BackwardEulerODESolver::assemble_rhs_vector(const MechanicalParams &    mechanical_parameters,
                                                 const MultiMatrixAccessor & matrix_accessor,
                                                 MultiVecDerivId & f_id,
                                                 SofaCaribou::Algebra::BaseVector *      f)
{

    const auto h = mechanical_parameters.dt();

    // 1. Clear the force and dforce vectors (F := 0, dF := 0)
    MechanicalResetForceVisitor(&mechanical_parameters, f_id,
                                false /* onlyMapped */)
    .execute(this->getContext());
    MechanicalResetForceVisitor(&mechanical_parameters, sofa::core::VecDerivId::dforce(),
                                false /* onlyMapped */)
            .execute(this->getContext());

    // 2. Go down in the current context tree calling `addForce` on every force field components,
    //    then go up from the leaves calling `applyJT` on every mechanical mappings
    //    NOTE: This correspond to the force terms dependant on the displacement, i.e. :
    //                             - f_0 = (-Ku + F)
    MechanicalComputeForceVisitor(&mechanical_parameters, f_id,
                                  true /* accumulate (to mapped node) */, true /*neglectingCompliance*/)
    .execute(this->getContext());

    // 3. Go down in the current context tree calling `addMBKdx` on every force field components,
    //    then go up from the leaves calling `applyJT` on every mechanical mappings. If the k
    //    mechanical parameter is not zero, the `applyDJT` will also be applied.
    //    NOTE: This correspond to the force terms dependant on the velocity, i.e. the damping terms :
    //                        - f_1 =  [- r_m M  - C  +  r_k K] v
    //          where K is in fact -K by SOFA's convention, hence the positive (+) sign.

    // Copy the mechanical parameters to temporarily swap the dx and v vectors
    // This is a hack since the BaseForcefield class doesn't take a multi-vector identified for the
    // dx in MBK[dx] multiplication. Since we want to do MBK[v], we need to swap the dv and the dx.
    auto m_params = mechanical_parameters;
    m_params.setDx(p_previous_v_id);

    // Set the Rayleigh's damping coefficients
    m_params.setMFactor(-d_rayleigh_mass.getValue());
    m_params.setKFactor(d_rayleigh_stiffness.getValue());

    // If another damping term is computed, it should be in the B matrix.
    m_params.setBFactor(-1);

    // Compute MBK [dv] and store the result in f_id
    MechanicalAddMBKdxVisitor(&m_params, f_id,
                              true /* Accumulate the contribution of mapped mo */ )
    .execute(this->getContext());

    // 4. Same thing as before, but this time computing
    //                        - f_2 =  [-(1 + h r_m) M  - h C  +  h r_k K ] a
    // where K is in fact -K by SOFA's convention, hence the positive (+) sign.
    m_params.setDx(p_a_id);

    // Set the Rayleigh's damping coefficients
    m_params.setMFactor(-(1 + h*d_rayleigh_mass.getValue()));
    m_params.setKFactor(h * d_rayleigh_stiffness.getValue());

    // If another damping term is computed, it should be in the B matrix.
    m_params.setBFactor(-h);

    // Compute MBK [dv] and store the result in f_id
    MechanicalAddMBKdxVisitor(&m_params, f_id,
                              true /* Accumulate the contribution of mapped mo */ )
            .execute(this->getContext());


    // 4. Calls the "projectResponse" method of every `BaseProjectiveConstraintSet` objects found in the current
    //    context tree. For example, the `FixedConstraints` component will set entries of fixed nodes to zero.
    MechanicalApplyConstraintsVisitor(&mechanical_parameters, f_id,
                                      nullptr /* W (also project the given compliance matrix) */)
    .execute(this->getContext());

    // 5. Copy force vectors from every top level (unmapped) mechanical objects into the given system vector f
    MechanicalMultiVectorToBaseVectorVisitor(&mechanical_parameters, f_id /* source */, f /* destination */, &matrix_accessor)
    .execute(this->getContext());
}

// Assemble A in A [da] = F
//          A = (1 + h*r_m) M   +   h C  +  [h * (h + r_k)] K
void BackwardEulerODESolver::assemble_system_matrix(const MechanicalParams & mechanical_parameters,
                                                    DefaultMultiMatrixAccessor & matrix_accessor,
                                                    SofaCaribou::Algebra::BaseMatrix * A)
{
    const auto h = mechanical_parameters.dt();

    // Step 1. Building stage:
    //         Here we go down on the current context sub-graph and call :
    //           1. ff->addKToMatrix(&K) and f->addBToMatrix() for every force field "ff" found.
    //           2. pc->applyConstraint(&K) for every BaseProjectiveConstraintSet "pc" found.
    //         If a mechanical mapping "m" is found during the traversal, and m->areMatricesMapped() is false, the
    //         traversal stops in the subgraph of the mapping.
    matrix_accessor.setGlobalMatrix(A);
    auto m_params = mechanical_parameters;
    m_params.setMFactor(1 + h*d_rayleigh_mass.getValue());
    m_params.setBFactor(h);
    m_params.setKFactor(-h*(h+d_rayleigh_stiffness.getValue())); // Here we multiply by -1 since K is in fact -K by SOFA's convention
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

// Propagate dv that was previously solved in A [dv] = F
void BackwardEulerODESolver::propagate_solution_increment(const MechanicalParams & mechanical_parameters,
                                                          const MultiMatrixAccessor & matrix_accessor,
                                                          const SofaCaribou::Algebra::BaseVector * dx,
                                                          MultiVecCoordId & x_id, MultiVecDerivId & v_id,
                                                          MultiVecDerivId & dx_id) {
    // 1. Prepare everything need to solve constraints later on
    const auto h = mechanical_parameters.dt();
    using Direction = sofa::core::objectmodel::BaseContext::SearchDirection;
    auto constraint_parameters = sofa::core::ConstraintParams(mechanical_parameters);
    auto constraint_solvers = this->getContext()->getObjects<sofa::core::behavior::ConstraintSolver>(Direction::Local);


    // 2. Copy vectors from the global system vector into every top level (unmapped) mechanical objects.
    MechanicalMultiVectorFromBaseVectorVisitor(&mechanical_parameters, dx_id, dx, &matrix_accessor).execute(this->getContext());

    // 3. a_{i+1}^n = a_{i}^n + da
    MechanicalVOpVisitor(&mechanical_parameters, p_a_id, p_a_id, dx_id).execute(this->getContext());

    // 4. v_{i+1}^n = v_0^n + h * a_{i+1}^n
    MechanicalVOpVisitor(&mechanical_parameters, v_id, p_previous_v_id, p_a_id, h).execute(this->getContext());

    // 5. Solve velocity constraints
    constraint_parameters.setOrder(sofa::core::ConstraintParams::VEL);
    for (auto * solver : constraint_solvers) {
        solver->solveConstraint(&constraint_parameters, v_id);
    }

    // 6. x_{i+1}^n = x_0^n + h v_{i+1}^n    with x_0^n the position vector at the beginning of the time step
    MechanicalVOpVisitor(&mechanical_parameters, x_id, p_previous_x_id, v_id, h).execute(this->getContext());

    // 7. Solve position constraints
    constraint_parameters.setOrder(sofa::core::ConstraintParams::POS);
    for (auto * solver : constraint_solvers) {
        solver->solveConstraint(&constraint_parameters, x_id);
    }

    // 8. Propagate positions to mapped mechanical objects, for example, identity mappings, barycentric mappings, etc.
    //    This will call the methods apply and applyJ on every mechanical mappings.
    MechanicalPropagateOnlyPositionAndVelocityVisitor(&mechanical_parameters).execute(this->getContext());

}

} // namespace SofaCaribou::ode
