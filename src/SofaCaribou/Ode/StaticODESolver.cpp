#include "StaticODESolver.h"

#include <sofa/core/ObjectFactory.h>
#include <sofa/helper/AdvancedTimer.h>
#include <sofa/simulation/MechanicalOperations.h>
#include <sofa/simulation/VectorOperations.h>
#include <sofa/simulation/PropagateEventVisitor.h>

namespace SofaCaribou {
namespace GraphComponents {
namespace ode {

using sofa::core::VecId;
using namespace sofa::defaulttype;
using namespace sofa::core::behavior;

StaticODESolver::StaticODESolver()
    : d_newton_iterations(initData(&d_newton_iterations,
            (unsigned) 1,
            "newton_iterations",
            "Number of newton iterations between each load increments (normally, one load increment per simulation time-step."))
    , d_correction_tolerance_threshold(initData(&d_correction_tolerance_threshold,
            (double) 1e-5,
            "correction_tolerance_threshold",
            "Convergence criterion: The newton iterations will stop when the norm of correction |du| reach this threshold."))
    , d_residual_tolerance_threshold( initData(&d_residual_tolerance_threshold,
            (double) 1e-5,
            "residual_tolerance_threshold",
            "Convergence criterion: The newton iterations will stop when the ratio between norm of the residual R_k = |f_k - K(u_k)| at iteration k over R_0 is lower than this threshold. "
            "Use a negative value to disable this criterion."))
    , d_shoud_diverge_when_residual_is_growing( initData(&d_shoud_diverge_when_residual_is_growing,
            false,
            "shoud_diverge_when_residual_is_growing",
            "Divergence criterion: The newton iterations will stop when the residual is greater than the one from the previous iteration."))
    , d_converged(initData(&d_converged, false, "converged", "Whether or not the last call to solve converged", true /*is_displayed_in_gui*/, true /*is_read_only*/))

{}

void StaticODESolver::solve(const sofa::core::ExecParams* params, double /*dt*/, sofa::core::MultiVecCoordId xResult, sofa::core::MultiVecDerivId /*vResult*/) {
    sofa::simulation::common::VectorOperations vop( params, this->getContext() );
    sofa::simulation::common::MechanicalOperations mop( params, this->getContext() );
    sofa::simulation::common::VisitorExecuteFunc executeVisitor(*this->getContext());

    MultiVecCoord x_start(&vop, sofa::core::VecCoordId::position() );
    MultiVecCoord x(&vop, xResult /*core::VecCoordId::position()*/ );
    MultiVecDeriv force( &vop, sofa::core::VecDerivId::force() );
    dx.realloc( &vop, true );
    U.realloc( &vop, true );

    // MO vector dx is not allocated by default, it will seg fault if the CG is used (dx is taken by default) with an IdentityMapping
    MultiVecDeriv tempdx(&vop, sofa::core::VecDerivId::dx() ); tempdx.realloc( &vop, true, true );

    // Set implicit param to true to trigger nonlinear stiffness matrix recomputation
    mop->setImplicit(true);

    msg_info() << "======= Starting static ODE solver in time step " << this->getTime();
    msg_info() << "(doing a maximum of " << d_newton_iterations.getValue() << " newton iterations)";

    unsigned n_it=0;
    double dx_norm = 0, du_norm = 0, R0 = 0, R = 0, Rn = 0;
    const auto & shoud_diverge_when_residual_is_growing = d_shoud_diverge_when_residual_is_growing.getValue();
    const auto & correction_tolerance_threshold = d_correction_tolerance_threshold.getValue();
    const auto & residual_tolerance_threshold = d_residual_tolerance_threshold.getValue();
    bool converged = false;
    const auto & newton_iterations = d_newton_iterations.getValue();

    sofa::helper::AdvancedTimer::stepBegin("StaticODESolver::Solve");

    while (n_it < newton_iterations) {
        sofa::helper::AdvancedTimer::stepBegin("NewtonStep");

        // If this is the first newton step, compute the first residual and make sure we didn't already reach the equilibrium.
        if (n_it == 0) {
            // compute addForce, in mapped: addForce + applyJT (vec)
            sofa::helper::AdvancedTimer::stepBegin("ComputeForce");

            // Reset the force vectors on every mechanical objects found in the current context tree
            // todo(jnbrunet): force.clear is probably not needed since mop.computeForce clears the forces by default
            force.clear();

            // Accumulate the force vectors
            // 1. Go down in the current context tree calling addForce on every forcefields
            // 2. Go up from the current context tree leaves calling applyJT on every mechanical mappings
            mop.computeForce(force);

            // Calls the "projectResponse" method of every BaseProjectiveConstraintSet objects found in the
            // current context tree.
            mop.projectResponse(force);
            sofa::helper::AdvancedTimer::stepEnd("ComputeForce");

            // Compute the initial residual
            R = sqrt(force.dot(force));
            R0 = R;

            if (residual_tolerance_threshold > 0 && R <= residual_tolerance_threshold) {
                msg_info() << "The ODE has already reached an equilibrium state";
                converged = true;
            }
        }

        sofa::helper::AdvancedTimer::stepBegin("MBKBuild");
        // 1. The MechanicalMatrix::K is an empty matrix that stores three floats called factors: m, b and k.
        // 2. the * operator simply multiplies each of the three factors a value. No matrix is built yet.
        // 3. The = operator first search for a linear solver in the current context. It then calls the "setSystemMBKMatrix"
        //    method of the linear solver.

        //    A. For LinearSolver using a GraphScatteredMatrix, nothing appends.
        //    B. For LinearSolver using other type of matrices (FullMatrix, SparseMatrix, CompressedRowSparseMatrix),
        //       the "addMBKToMatrix" method is called on each BaseForceField objects and the "applyConstraint" method
        //       is called on every BaseProjectiveConstraintSet objects.
        sofa::core::behavior::MultiMatrix<sofa::simulation::common::MechanicalOperations> matrix(&mop);
        matrix = MechanicalMatrix::K * -1.0;
        sofa::helper::AdvancedTimer::stepEnd("MBKBuild");

        if (not converged) {
            // Solving the system
            // for CG: calls iteratively addDForce, mapped:  [applyJ, addDForce, applyJt(vec)]+
            // for LDL: solves the system, everything's already assembled
            sofa::helper::AdvancedTimer::stepBegin("MBKSolve");

            // Calls methods "setSystemRHVector", "setSystemLHVector" and "solveSystem" of the LinearSolver component
            matrix.solve(dx, force);
            sofa::helper::AdvancedTimer::stepEnd("MBKSolve");

            // Updating the geometry
            x.eq(x_start, dx, 1);

            // Solving constraints
            // Calls "solveConstraint" method of every ConstraintSolver objects found in the current context tree.
            mop.solveConstraint(x, sofa::core::ConstraintParams::POS);

            //propagate positions to mapped nodes (calls apply, applyJ)
            sofa::core::MechanicalParams mp;
            sofa::simulation::MechanicalPropagateOnlyPositionAndVelocityVisitor(&mp).execute(this->getContext());

            // compute addForce, in mapped: addForce + applyJT (vec)
            sofa::helper::AdvancedTimer::stepBegin("ComputeForce");
            force.clear();
            mop.computeForce(force);
            mop.projectResponse(force);
            sofa::helper::AdvancedTimer::stepEnd("ComputeForce");

            // Residual
            R = sqrt(force.dot(force));

            if (n_it == 0) {
                R0 = R;
            }

            // Displacement
            U.peq(dx);
            dx_norm = sqrt(dx.dot(dx));
            du_norm = sqrt(U.dot(U));
        }

        if (not converged) {
            msg_info() << "Newton iteration #" << n_it + 1
                       << ": |R|/|R0| = " << R/R0    << " (threshold of " << residual_tolerance_threshold   << ")"
                       << "  |du| = "     << dx_norm << " (threshold of " << correction_tolerance_threshold << ")";
        }
        sofa::helper::AdvancedTimer::valSet("residual", R);
        sofa::helper::AdvancedTimer::valSet("correction", dx_norm);
        sofa::helper::AdvancedTimer::valSet("displacement", du_norm);
        sofa::helper::AdvancedTimer::stepEnd("NewtonStep");

        if (not converged and shoud_diverge_when_residual_is_growing and R > Rn and n_it > 1) {
            msg_info() << "[DIVERGED] Residual's norm increased from "<< Rn << " to " << R;
        }

        if (not converged and correction_tolerance_threshold > 0 and dx_norm < correction_tolerance_threshold*du_norm) {
            converged = true;
            msg_info() << "[CONVERGED] The correction's ratio |du|/|U| = " << dx_norm/du_norm << " is smaller than the threshold of "
                       << correction_tolerance_threshold;
        }

        if (not converged and residual_tolerance_threshold > 0 and R < residual_tolerance_threshold*R0) {
            converged = true;
            msg_info() << "[CONVERGED] The residual's ratio |R|/|R0| = " << R/R0 << " is smaller than the threshold of "
                       << residual_tolerance_threshold;
        }


        // Save the last residual to check the growing residual criterion at the next step
        Rn = R;

        if (converged) {
            break;
        }

        n_it++;

    } // End while (n_it < d_newton_iterations.getValue())

    if (n_it >= newton_iterations) {
        n_it--;
        msg_info() << "[DIVERGED] The number of Newton iterations reached the maximum of " << newton_iterations << " iterations";
    }

    d_converged.setValue(converged);

    sofa::helper::AdvancedTimer::valSet("has_converged", converged ? 1 : 0);
    sofa::helper::AdvancedTimer::valSet("nb_iterations", n_it+1);
    sofa::helper::AdvancedTimer::valSet("residual", Rn);
    sofa::helper::AdvancedTimer::valSet("correction", dx_norm);
    sofa::helper::AdvancedTimer::valSet("displacement", du_norm);
    sofa::helper::AdvancedTimer::stepEnd("StaticODESolver::Solve");
}


SOFA_DECL_CLASS(StaticODESolver)

int StaticODESolverClass = sofa::core::RegisterObject("Static ODE Solver")
    .add< StaticODESolver >()
;

} // namespace ode
} // namespace GraphComponents
} // namespace SofaCaribou
