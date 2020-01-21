#include "ConjugateGradientSolver.h"
#include <sofa/core/ObjectFactory.h>
#include <sofa/helper/AdvancedTimer.h>
#include <sofa/simulation/MechanicalOperations.h>
#include <sofa/simulation/VectorOperations.h>

namespace SofaCaribou::GraphComponents::solver {

ConjugateGradientSolver::ConjugateGradientSolver()
: d_maximum_number_of_iterations(initData(&d_maximum_number_of_iterations,
    (unsigned int) 25,
    "maximum_number_of_iterations",
    "Maximum number of iterations before diverging."))
, d_residual_tolerance_threshold(initData(&d_residual_tolerance_threshold,
    1e-5,
    "residual_tolerance_threshold",
    "Convergence criterion: The CG iterations will stop when the ratio between norm of the residual "
    "r(k+1) = |r(k) - a(k) A p(k) | at iteration k+1 over r_0 is lower than this threshold."))
{
}

void ConjugateGradientSolver::setSystemMBKMatrix(const sofa::core::MechanicalParams* mparams) {
    p_mechanical_params = mparams;
}

void ConjugateGradientSolver::setSystemRHVector(sofa::core::MultiVecDerivId b_id) {
    p_b_id = b_id;
}

void ConjugateGradientSolver::setSystemLHVector(sofa::core::MultiVecDerivId x_id) {
    p_x_id = x_id;
}

void ConjugateGradientSolver::solveSystem() {
    sofa::simulation::common::VectorOperations vop( p_mechanical_params, this->getContext() );
    sofa::simulation::common::MechanicalOperations mop( p_mechanical_params, this->getContext() );

    // Create temporary vectors needed for the method
    MultiVecDeriv r(&vop);
    MultiVecDeriv p(&vop);
    MultiVecDeriv q(&vop); // q = A*p

    // Gather the x and b vector identifiers
    MultiVecDeriv x(&vop, p_x_id);
    MultiVecDeriv b(&vop, p_b_id);

    // Get the method parameters
    const auto & maximum_number_of_iterations = d_maximum_number_of_iterations.getValue();
    const auto & residual_tolerance_threshold = d_residual_tolerance_threshold.getValue();

    // Get the matrices coefficient m, b and k : A = (mM + bB + kK)
    const auto  m_coef = p_mechanical_params->mFactor();
    const auto  b_coef = p_mechanical_params->bFactor();
    const auto  k_coef = p_mechanical_params->kFactor();

    // Declare the method variables
    FLOATING_POINT_TYPE b_norm, r_norm; // Residual norm
    FLOATING_POINT_TYPE rho0, rho1; // Stores r*r as it is used two times per iterations
    FLOATING_POINT_TYPE alpha, beta; // Alpha and Beta coefficients
    UNSIGNED_INTEGER_TYPE iteration_number = 0; // Current iteration number

    sofa::helper::AdvancedTimer::stepBegin("ConjugateGradient::solve");

    // Compute the first residual r(0) = b - A*x(0) = b - q
    mop.propagateDxAndResetDf(x, q); // Set q = 0 and calls applyJ(x) on every mechanical mappings
    mop.addMBKdx(q, m_coef, b_coef, k_coef, false); // q = (m M + b B + k K) x
    mop.projectResponse(q); // BaseProjectiveConstraintSet::projectResponse(q)
    r.eq( b, q, -1.0 );   // r = b - q

    // Check for initial convergence: |r0|/|b| < threshold
    rho0 = r.dot(r);
    r_norm = sqrt(rho0);
    b_norm = b.norm();
    if (r_norm < residual_tolerance_threshold*b_norm) {
        msg_info() << "The linear system has already reached an equilibrium state";
        msg_info() << "|R| = " << r_norm << ", |b| = " << b_norm << ", threshold = " << residual_tolerance_threshold;
    } else {
        p = r; // p(0) = r(0)
        while (iteration_number < maximum_number_of_iterations) {
            // 1. Computes q(k+1) = A*p(k)
            mop.propagateDxAndResetDf(p, q); // Set q = 0 and calls applyJ(p) on every mechanical mappings
            mop.addMBKdx(q, m_coef, b_coef, k_coef, false); // q = (m M + b B + k K) x
            mop.projectResponse(q); // BaseProjectiveConstraintSet::projectResponse(q)

            // 2. Computes x(k+1) and r(k+
            alpha = rho0 / p.dot(q);
            x.peq(p,alpha);  // x = x + alpha*p
            r.peq(q,-alpha); // r = r - alpha*q

            // 3. Computes the new residual norm
            rho1 = r.dot(r);
            r_norm = sqrt(rho1);

            // 4. Print information on the current iteration
            msg_info() << "CG iteration #" << iteration_number+1
                       << ": |r0|/|b| = "  << r_norm/b_norm
                       << ", threshold = " << residual_tolerance_threshold;

            // 5. Check for convergence: |r|/|b| < threshold
            if (r_norm< residual_tolerance_threshold*b_norm) {
                msg_info() << "CG converged!";
                break;
            }

            // 6. Compute p(k+1)
            beta = rho1 / rho0;
            p.eq(r,p,beta); // p = r + beta*p

            rho0 = rho1;
            ++iteration_number;
        }
    }

    sofa::helper::AdvancedTimer::stepEnd("ConjugateGradient::solve");
}

int CGLinearSolverClass = sofa::core::RegisterObject("Linear system solver using the conjugate gradient iterative algorithm")
                              .add< ConjugateGradientSolver>(true);

} // namespace SofaCaribou::GraphComponents::solver