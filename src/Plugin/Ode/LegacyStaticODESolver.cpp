#include "LegacyStaticODESolver.h"

#include <iomanip>
#include <chrono>

#include <SofaCaribou/Solver/ConjugateGradientSolver.h>

#include <sofa/core/ObjectFactory.h>
#include <sofa/helper/AdvancedTimer.h>
#include <sofa/simulation/MechanicalOperations.h>
#include <sofa/simulation/VectorOperations.h>
#include <sofa/simulation/Node.h>

namespace SofaCaribou::ode {

using sofa::core::VecId;
using namespace sofa::defaulttype;
using namespace sofa::core::behavior;

LegacyStaticODESolver::LegacyStaticODESolver()
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
    , d_warm_start(initData(&d_warm_start,
        false,
        "warm_start",
        "For iterative linear solvers, use the previous solution has a warm start. "
        "Note that for the first newton step, the current position is used as the warm start."))
    , d_converged(initData(&d_converged, false, "converged", "Whether or not the last call to solve converged", true /*is_displayed_in_gui*/, true /*is_read_only*/))
{}

void LegacyStaticODESolver::solve(const sofa::core::ExecParams* params, double /*dt*/, sofa::core::MultiVecCoordId xResult, sofa::core::MultiVecDerivId /*vResult*/) {
    using namespace sofa::helper::logging;
    using namespace std::chrono;
    using std::chrono::steady_clock;

    sofa::core::MechanicalParams mparams (*params);
    mparams.setDx(dx);

    sofa::simulation::common::VectorOperations vop( &mparams, this->getContext() );
    sofa::simulation::common::MechanicalOperations mop( &mparams, this->getContext() );

    MultiVecCoord x(&vop, xResult);
    MultiVecDeriv force( &vop, sofa::core::VecDerivId::force() );

    // Options
    const auto & shoud_diverge_when_residual_is_growing = d_shoud_diverge_when_residual_is_growing.getValue();
    const auto & correction_tolerance_threshold = d_correction_tolerance_threshold.getValue();
    const auto & residual_tolerance_threshold = d_residual_tolerance_threshold.getValue();
    const auto & newton_iterations = d_newton_iterations.getValue();
    const auto & warm_start = d_warm_start.getValue();
    const auto & print_log = f_printLog.getValue();
    auto info = MessageDispatcher::info(Message::Runtime, ComponentInfo::SPtr(new ComponentInfo(this->getClassName())), SOFA_FILE_INFO);

    // Get the linear solver and convert it to a CG if it is one
    const auto context = mop.ctx;
    auto linear_solver = context->get<LinearSolver>(context->getTags(), sofa::core::objectmodel::BaseContext::SearchDown);
    auto cg_linear_solver = dynamic_cast<SofaCaribou::solver::ConjugateGradientSolver *>(linear_solver);

    // Incremental displacement of one iteration
    dx.realloc( &vop );

    if (not warm_start) {
        // dx := 0
        dx.clear();
    }

    // Total displacement increment since the beginning
    U.realloc( &vop );
    U.clear();

    // Set implicit param to true to trigger nonlinear stiffness matrix recomputation
    mop->setImplicit(true);

    if (print_log) {
        info << "======= Starting static ODE solver =======\n";
        info << "Time step             : " << this->getTime() << "\n";
        info << "Context               : " << dynamic_cast<sofa::simulation::Node *>(this->getContext())->getPathName() << "\n";
        info << "Max iterations        : " << newton_iterations << "\n";
        info << "Residual tolerance    : " << residual_tolerance_threshold << "\n";
        info << "Correction tolerance  : " << correction_tolerance_threshold << "\n\n";
    }

    unsigned n_it=0;
    double dx_squared_norm, du_squared_norm, R_squared_norm = 0, Rn_squared_norm = 0;
    const auto squared_residual_threshold = residual_tolerance_threshold*residual_tolerance_threshold;
    const auto squared_correction_threshold = correction_tolerance_threshold*correction_tolerance_threshold;
    bool converged = false, diverged = false;
    steady_clock::time_point t;

    // Resize vectors containing the newton residual norms
    p_squared_residuals.clear();
    p_squared_residuals.reserve(newton_iterations);

    // Resize vectors containing the times took to compute the newton iterations
    p_times.clear();
    p_times.reserve(newton_iterations);

    // If the linear solver is a CG, resize the vector containing CG residual norms
    if (cg_linear_solver) {
        p_iterative_linear_solver_squared_residuals.clear();
        p_iterative_linear_solver_squared_rhs_norms.clear();
        p_iterative_linear_solver_squared_residuals.reserve(newton_iterations);
        p_iterative_linear_solver_squared_rhs_norms.reserve(newton_iterations);
    }


    sofa::helper::AdvancedTimer::stepBegin("StaticODESolver::Solve");

    // #########################
    //      First residual
    // #########################
    // Before starting the newton iterations, we first need to compute
    // the residual with the updated right-hand side (the new load increment)

    // compute addForce, in mapped: addForce + applyJT (vec)
    sofa::helper::AdvancedTimer::stepBegin("ComputeForce");

    // Accumulate the force vectors
    // 1. Clear the force vector (F := 0)
    // 2. Go down in the current context tree calling addForce on every forcefields
    // 3. Go up from the current context tree leaves calling applyJT on every mechanical mappings
    mop.computeForce(force);

    // Calls the "projectResponse" method of every BaseProjectiveConstraintSet objects found in the
    // current context tree. An example of such constraint set is the FixedConstraint. In this case,
    // it will set to 0 every row (i, _) of the right-hand side (force) vector for the ith degree of
    // freedom.
    mop.projectResponse(force);
    sofa::helper::AdvancedTimer::stepEnd("ComputeForce");

    // Compute the initial residual
    R_squared_norm = force.dot(force);
    p_squared_initial_residual = R_squared_norm;

    if (residual_tolerance_threshold > 0 && R_squared_norm <= residual_tolerance_threshold) {
        converged = true;
        if (print_log) {
            info << "The ODE has already reached an equilibrium state" << "\n";
        }
    }

    // #########################
    // Newton-Raphson iterations
    // #########################
    while (not converged and n_it < newton_iterations) {
        sofa::helper::AdvancedTimer::stepBegin("NewtonStep");

        t = steady_clock::now();

        sofa::helper::AdvancedTimer::stepBegin("MBKBuild");
        // 1. The MechanicalMatrix::K is an empty matrix that stores three floats called factors: m, b and k.
        // 2. the * operator simply multiplies each of the three factors a value. No matrix is built yet.
        // 3. The = operator first search for a linear solver in the current context. It then calls the "setSystemMBKMatrix"
        //    method of the linear solver.

        //    A. For LinearSolver using a GraphScatteredMatrix (ie, non-assembled matrices), nothing appends.
        //    B. For LinearSolver using other type of matrices (FullMatrix, SparseMatrix, CompressedRowSparseMatrix),
        //       the "addMBKToMatrix" method is called on each BaseForceField objects and the "applyConstraint" method
        //       is called on every BaseProjectiveConstraintSet objects. An example of such constraint set is the
        //       FixedConstraint. In this case, it will set to 0 every column (_, i) and row (i, _) of the assembled
        //       matrix for the ith degree of freedom.
        sofa::core::behavior::MultiMatrix<sofa::simulation::common::MechanicalOperations> matrix(&mop);
        matrix = MechanicalMatrix::K * -1.0;
        sofa::helper::AdvancedTimer::stepEnd("MBKBuild");

        // Solving the system
        // for CG: calls iteratively addDForce, mapped:  [applyJ, addDForce, applyJt(vec)]+
        // for LDL: solves the system, everything's already assembled
        sofa::helper::AdvancedTimer::stepBegin("MBKSolve");

        // Calls methods "setSystemRHVector", "setSystemLHVector" and "solveSystem" of the LinearSolver component
        matrix.solve(dx, force);
        sofa::helper::AdvancedTimer::stepEnd("MBKSolve");

        // Updating the geometry
        x.peq(dx); // x := x + dx

        // Solving constraints
        // Calls "solveConstraint" method of every ConstraintSolver objects found in the current context tree.
        // todo(jnbrunet): Shouldn't this be done AFTER the position propagation of the mapped nodes?
        mop.solveConstraint(x, sofa::core::ConstraintParams::POS);

        // Propagate positions to mapped mechanical objects, for example, identity mappings, barycentric mappings, ...
        // This will call the methods apply and applyJ on every mechanical mappings.
        sofa::simulation::MechanicalPropagateOnlyPositionAndVelocityVisitor(&mparams).execute(this->getContext());

        // Update the force residual
        sofa::helper::AdvancedTimer::stepBegin("ComputeForce");

        // Accumulate the force vectors
        // 1. Clear the force vector (F := 0)
        // 2. Go down in the current context tree calling addForce on every forcefields
        // 3. Go up from the current context tree leaves calling applyJT on every mechanical mappings
        mop.computeForce(force);

        // Calls the "projectResponse" method of every BaseProjectiveConstraintSet objects found in the
        // current context tree.
        mop.projectResponse(force);
        sofa::helper::AdvancedTimer::stepEnd("ComputeForce");

        // Residual
        R_squared_norm = force.dot(force);

        // Displacement
        U.peq(dx);
        dx_squared_norm = dx.dot(dx);
        du_squared_norm= U.dot(U);

        // Stop the timers here as the computation is done
        auto iteration_time = duration_cast<nanoseconds>(steady_clock::now() - t).count();
        p_times.emplace_back(static_cast<UNSIGNED_INTEGER_TYPE>(iteration_time));
        sofa::helper::AdvancedTimer::stepEnd("NewtonStep");

        p_squared_residuals.emplace_back(R_squared_norm);

        // If the linear solver is the caribou's CG, copy the residual norms of the CG iterations
        if (cg_linear_solver) {
            p_iterative_linear_solver_squared_residuals.emplace_back(cg_linear_solver->squared_residuals());
            p_iterative_linear_solver_squared_rhs_norms.emplace_back(cg_linear_solver->squared_initial_residual());
        }

        // We completed one iteration, increment the counter
        n_it++;

        if( print_log ) {
            info << "Newton iteration #" << std::left << std::setw(5)  << n_it
                 << std::scientific
                 << "  |R|/|R0| = " << std::setw(12) << sqrt(R_squared_norm  / p_squared_residuals[0])
                 << "  |du| = "     << std::setw(12) << sqrt(dx_squared_norm / du_squared_norm)
                 << std::defaultfloat;
            if (cg_linear_solver) {
                info << "  CG iterations = " << std::setw(5) << cg_linear_solver->squared_residuals().size();
            }
            info << "  Time = " << iteration_time/1000/1000 << " ms";
            info << "\n";
        }

        if (std::isnan(R_squared_norm) or std::isnan(dx_squared_norm) or du_squared_norm < EPSILON) {
            diverged = true;
            if (print_log) {
                info << "[DIVERGED]";
                if (std::isnan(R_squared_norm)) {
                    info << " The residual's ratio |R| is NaN.";
                }
                if (std::isnan(dx_squared_norm)) {
                    info << " The correction's ratio |du| is NaN.";
                }
                if (du_squared_norm < EPSILON) {
                    info << " The correction's ratio |du|/|U| is NaN (|U| is zero).";
                }
                info << "\n";
            }
            break;
        }


        if (correction_tolerance_threshold > 0 and dx_squared_norm < squared_correction_threshold*du_squared_norm) {
            converged = true;
            if (print_log) {
                info  << "[CONVERGED] The correction's ratio |du|/|U| = " << sqrt(dx_squared_norm/du_squared_norm) << " is smaller than the threshold of " << correction_tolerance_threshold << ".\n";
            }
            break;
        }

        if (residual_tolerance_threshold > 0 and R_squared_norm < squared_residual_threshold*p_squared_residuals[0]) {
            converged = true;
            if (print_log) {
                info << "[CONVERGED] The residual's ratio |R|/|R0| = " << sqrt(R_squared_norm/p_squared_residuals[0]) << " is smaller than the threshold of " << residual_tolerance_threshold << ".\n";
            }
            break;
        }

        if (shoud_diverge_when_residual_is_growing and R_squared_norm > Rn_squared_norm and n_it > 1) {
            diverged = true;
            if (print_log) {
                info << "[DIVERGED] Residual's norm increased from "<< sqrt(Rn_squared_norm) << " to " << sqrt(R_squared_norm) << ".\n";
            }
            break;
        }

        // Save the last residual to check the growing residual criterion at the next step
        Rn_squared_norm = R_squared_norm;

        if (not warm_start) {
            dx.clear();
        }

    } // End while (not converged and not diverged and n_it < newton_iterations)

    n_it--; // Reset to the actual index of the last iteration completed

    if (not converged and not diverged and n_it == (newton_iterations-1)) {
        if (print_log) {
            info << "[DIVERGED] The number of Newton iterations reached the maximum of " << newton_iterations << " iterations" << ".\n";
        }
    }

    d_converged.setValue(converged);

    sofa::helper::AdvancedTimer::valSet("has_converged", converged ? 1 : 0);
    sofa::helper::AdvancedTimer::valSet("nb_iterations", n_it+1);
    sofa::helper::AdvancedTimer::stepEnd("StaticODESolver::Solve");
}


SOFA_DECL_CLASS(StaticODESolver)

int StaticODESolverClass = sofa::core::RegisterObject("Static ODE Solver")
    .add< LegacyStaticODESolver >()
;

} // namespace SofaCaribou::ode
