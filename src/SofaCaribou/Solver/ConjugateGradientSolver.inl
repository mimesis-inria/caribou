#include <SofaCaribou/Solver/ConjugateGradientSolver.h>
#include <SofaCaribou/Solver/EigenSolver.inl>
#include<SofaCaribou/Algebra/EigenMatrix.h>
#include <SofaCaribou/Visitor/AssembleGlobalMatrix.h>
#include <SofaCaribou/Visitor/ConstrainGlobalMatrix.h>
#include <Caribou/macros.h>

DISABLE_ALL_WARNINGS_BEGIN
#include <sofa/helper/AdvancedTimer.h>
#if (defined(SOFA_VERSION) && SOFA_VERSION >= 210600)
#include <sofa/helper/ScopedAdvancedTimer.h>
#endif
#include <sofa/simulation/MechanicalOperations.h>
#include <sofa/simulation/VectorOperations.h>
DISABLE_ALL_WARNINGS_END

#include <iomanip>

#if !EIGEN_VERSION_AT_LEAST(3,3,0)
namespace Eigen {
using Index = EIGEN_DEFAULT_DENSE_INDEX_TYPE;
}
#endif

namespace SofaCaribou::solver {

using Timer = sofa::helper::AdvancedTimer;
using Algebra::EigenMatrix;

template <class EigenMatrix_t>
ConjugateGradientSolver<EigenMatrix_t>::ConjugateGradientSolver()
: d_verbose(initData(&d_verbose,
    false,
    "verbose",
    "When the data attribute 'printLog' is activated, activating the verbose flag will make the CG output convergence "
    "information at every iterations in addition to the usual information printed."))
, d_maximum_number_of_iterations(initData(&d_maximum_number_of_iterations,
    (unsigned int) 25,
    "maximum_number_of_iterations",
    "Maximum number of iterations before diverging."))
, d_residual_tolerance_threshold(initData(&d_residual_tolerance_threshold,
    1e-5,
    "residual_tolerance_threshold",
    "Convergence criterion: The CG iterations will stop when the ratio between norm of the residual "
    "r(k+1) = |r(k) - a(k) A p(k) | at iteration k+1 over r_0 is lower than this threshold."))
, d_preconditioning_method(initData(&d_preconditioning_method,
    "preconditioning_method",
    R"(
        Preconditioning methods are:
            None:                No preconditioning, hence the complete matrix won't be built. (default)
            Identity:            A naive preconditioner which approximates any matrix as the identity matrix.
            Diagonal:            Preconditioning using an approximation of A.x = b by ignoring all off-diagonal entries of A.
    )"
#if EIGEN_VERSION_AT_LEAST(3,3,0)
    R"(
            IncompleteCholesky:  Preconditioning based on the incomplete Cholesky factorization.
    )"
#endif
    R"(
            IncompleteLU:        Preconditioning based on the incomplete LU factorization.
    )",
    true /*displayed_in_GUI*/, false /*read_only_in_GUI*/))
{
    // Explicitly state the available preconditioning methods
    p_preconditioners.emplace_back("None", PreconditioningMethod::None);
    p_preconditioners.emplace_back("Identity", PreconditioningMethod::Identity);
    p_preconditioners.emplace_back("Diagonal", PreconditioningMethod::Diagonal);
#if EIGEN_VERSION_AT_LEAST(3,3,0)
    p_preconditioners.emplace_back("IncompleteCholesky", PreconditioningMethod::IncompleteCholesky);
#endif
    p_preconditioners.emplace_back("IncompleteLU", PreconditioningMethod::IncompleteLU);

    // Fill-in the data option group with the available preconditioning methods
    std::vector<std::string> preconditioner_names;
    for (const auto & preconditioner : p_preconditioners) {
        const std::string & preconditioner_name = preconditioner.first;
        preconditioner_names.emplace_back(preconditioner_name);
    }
    d_preconditioning_method.setValue(sofa::helper::OptionsGroup(preconditioner_names));
    sofa::helper::WriteAccessor<Data< sofa::helper::OptionsGroup >> preconditioning_method = d_preconditioning_method;
    preconditioning_method->setSelectedItem((unsigned int) 1);
}

template <class EigenMatrix_t>
auto ConjugateGradientSolver<EigenMatrix_t>::get_preconditioning_method_from_string(const std::string & preconditioner_name) const -> PreconditioningMethod{
    auto to_lower = [](const std::string & input_string) {
        std::string output_string = input_string;
        std::transform(input_string.begin(), input_string.end(), output_string.begin(),
            [](unsigned char c){ return std::tolower(c); }
        );
        return output_string;
    };

    const std::string query = to_lower(preconditioner_name);

    for (const auto & preconditioner : p_preconditioners) {
        const std::string & preconditioner_name = preconditioner.first;
        const PreconditioningMethod & preconditioner_id = preconditioner.second;
        if (to_lower(preconditioner_name) == query) {
            return preconditioner_id;
        }
    }

    return PreconditioningMethod::None;
}

template <class EigenMatrix_t>
void ConjugateGradientSolver<EigenMatrix_t>::setSystemMBKMatrix(const sofa::core::MechanicalParams* mparams) {
    // Save the current mechanical parameters (m, b and k factors of the mass (M), damping (B) and
    // stiffness (K) matrices)
    p_mechanical_params = *mparams;

    // Get the preconditioning method
    const PreconditioningMethod preconditioning_method = get_preconditioning_method_from_string(d_preconditioning_method.getValue().getSelectedItem());

    // If we have a preconditioning method, the global system matrix has to be constructed from the current context subgraph.
    if (preconditioning_method != PreconditioningMethod::None) {
        Base::setSystemMBKMatrix(mparams);
    }
}

template <class EigenMatrix_t>
void ConjugateGradientSolver<EigenMatrix_t>::setSystemRHVector(sofa::core::MultiVecDerivId b_id) {
    // Register the RHS to the mechanical parameters
    p_mechanical_params.setDf(b_id);

    sofa::simulation::common::MechanicalOperations mop(&p_mechanical_params, this->getContext());
    p_b_id = b_id;

    // Get the preconditioning method
    const PreconditioningMethod preconditioning_method = get_preconditioning_method_from_string(d_preconditioning_method.getValue().getSelectedItem());

    // If we have a preconditioning method, we copy the vectors of the mechanical objects into a global eigen vector.
    if (preconditioning_method != PreconditioningMethod::None) {
        Base::setSystemRHVector(b_id);
    }
}

template <class EigenMatrix_t>
void ConjugateGradientSolver<EigenMatrix_t>::setSystemLHVector(sofa::core::MultiVecDerivId x_id) {
    // Register the dx vector to the mechanical parameters
    p_mechanical_params.setDx(x_id);

    sofa::simulation::common::MechanicalOperations mop(&p_mechanical_params, this->getContext());
    p_x_id = x_id;

    // Get the preconditioning method
    const PreconditioningMethod preconditioning_method = get_preconditioning_method_from_string(d_preconditioning_method.getValue().getSelectedItem());

    // If we have a preconditioning method, we copy the vectors of the mechanical objects into a global eigen vector.
    if (preconditioning_method != PreconditioningMethod::None) {
        Base::setSystemLHVector(x_id);
    }
}

template <class EigenMatrix_t>
void ConjugateGradientSolver<EigenMatrix_t>::solveSystem() {
    // Get the preconditioning method
    const PreconditioningMethod preconditioning_method = get_preconditioning_method_from_string(d_preconditioning_method.getValue().getSelectedItem());

    if (preconditioning_method == PreconditioningMethod::None) {
        sofa::helper::ScopedAdvancedTimer _t_("ConjugateGradient::solve");
        // Gather the x and b vector identifiers
        sofa::simulation::common::VectorOperations vop( &p_mechanical_params, this->getContext() );
        MultiVecDeriv x(&vop, p_x_id);
        MultiVecDeriv b(&vop, p_b_id);

        // Solve without having assembled the global matrix A (not needed since no preconditioning)
        solve(b, x);
    } else {
        // Solve using a preconditioning method. Here the global matrix A and the vectors x and b have been assembled.
        Base::solveSystem();
    }

}

template <class EigenMatrix_t>
void ConjugateGradientSolver<EigenMatrix_t>::solve(sofa::core::behavior::MultiVecDeriv & b, sofa::core::behavior::MultiVecDeriv & x) {
    sofa::simulation::common::VectorOperations vop( &p_mechanical_params, this->getContext() );
    sofa::simulation::common::MechanicalOperations mop( &p_mechanical_params, this->getContext() );

    // Create temporary vectors needed for the method
    sofa::core::behavior::MultiVecDeriv r(&vop);
    sofa::core::behavior::MultiVecDeriv p(&vop);
    sofa::core::behavior::MultiVecDeriv q(&vop);

    // Get the method parameters
    const auto & maximum_number_of_iterations = d_maximum_number_of_iterations.getValue();
    const auto & residual_tolerance_threshold = d_residual_tolerance_threshold.getValue();
    const auto & verbose = d_verbose.getValue();

    p_squared_residuals.clear();
    p_squared_residuals.reserve(maximum_number_of_iterations);

    // Get the matrices coefficient m, b and k : A = (mM + bB + kK)
    const auto  m_coef = p_mechanical_params.mFactor();
    const auto  b_coef = p_mechanical_params.bFactor();
    const auto  k_coef = p_mechanical_params.kFactor();

    // Declare the method variables
    FLOATING_POINT_TYPE b_norm_2 = 0., r_norm_2 = 0.; // RHS and residual squared norms
    FLOATING_POINT_TYPE rho0, rho1 = 0.; // Stores r*r as it is used two times per iterations
    FLOATING_POINT_TYPE alpha, beta; // Alpha and Beta coefficients
    FLOATING_POINT_TYPE threshold; // Residual threshold
    UNSIGNED_INTEGER_TYPE iteration_number = 0; // Current iteration number
    bool converged = false;
    const auto zero = (std::numeric_limits<FLOATING_POINT_TYPE>::min)(); // A numerical floating point zero

    // Make sure that the right hand side isn't zero
    b_norm_2 = b.dot(b);
    p_squared_initial_residual = b_norm_2;
    if (b_norm_2 < EPSILON) {
        msg_info() << "Right-hand side of the system is zero, hence x = 0.";
        x.clear();
        goto end; // The goto is important to catch the last timer call before ending the function
    }

    // Compute the tolerance w.r.t |b| since |r|/|b| < threshold is equivalent to  r^2 < b^2 * threshold^2
    // threshold = b^2 * residual_tolerance_threshold^2
    threshold = std::max(residual_tolerance_threshold*residual_tolerance_threshold*b_norm_2, zero);

    // INITIAL RESIDUAL
    // Do the A*x(0) with visitors since we did not construct the matrix A
    mop.propagateDxAndResetDf(x, q); // Set q = 0 and calls applyJ(x) on every mechanical mappings
    mop.addMBKdx(q, m_coef, b_coef, k_coef, false); // q = (m M + b B + k K) x

    // Do the projection of the result in the constrained space
    mop.projectResponse(q); // BaseProjectiveConstraintSet::projectResponse(q)

    // Finally, compute the initial residual r = b - A*x
    r.eq( b, q, -1.0 );   // r = b - q

    // Check for initial convergence: |r0|/|b| < threshold
    r_norm_2 = r.dot(r);
    if (r_norm_2 < threshold) {
        msg_info() << "The linear system has already reached an equilibrium state";
        msg_info() << "|r|/|b| = " << sqrt(r_norm_2/b_norm_2) << ", threshold = " << residual_tolerance_threshold;
        goto end; // The goto is important to catch the last timer call before ending the function
    }

    // Compute the initial search direction
    rho0 = r_norm_2;
    p = r; // p(0) = r(0)

    // ITERATIONS
    while (not converged and iteration_number < maximum_number_of_iterations) {
        Timer::stepBegin("cg_iteration");
        // 1. Computes q(k+1) = A*p(k)
        mop.propagateDxAndResetDf(p, q); // Set q = 0 and calls applyJ(p) on every mechanical mappings
        mop.addMBKdx(q, m_coef, b_coef, k_coef, false); // q = (m M + b B + k K) x

        // We need to project the residual in the constrained space since the constraints haven't been added to the
        // matrix (the matrix is never constructed) and addMBKdx of the forcefields do not take constraints into account.
        mop.projectResponse(q); // BaseProjectiveConstraintSet::projectResponse(q)

        // 2. Computes x(k+1) and r(k+1)
        alpha = rho0 / p.dot(q);
        x.peq(p,alpha);  // x = x + alpha*p
        r.peq(q,-alpha); // r = r - alpha*q

        // 3. Computes the new residual norm
        r_norm_2 = r.dot(r);
        p_squared_residuals.emplace_back(r_norm_2);

        // 4. Print information on the current iteration
        msg_info_when(verbose) << "CG iteration #" << iteration_number+1
                               << ": |r|/|b| = "   << sqrt(r_norm_2/b_norm_2)
                               << "(threshold is " << residual_tolerance_threshold << ")";

        // 5. Check for convergence: |r|/|b| < threshold
        if (r_norm_2 < threshold) {
            converged = true;
        } else {
            // 6. Compute the next search direction
            rho1 = r_norm_2;
            beta = rho1 / rho0;
            p.eq(r,p,beta); // p = r + beta*p

            rho0 = rho1;
        }

        ++iteration_number;
        Timer::stepEnd("cg_iteration");
    }

    iteration_number--; // Reset to the actual index of the last iteration completed

    if (converged) {
        msg_info() << "CG converged in " << (iteration_number+1)
                   << " iterations with a residual of |r|/|b| = " << sqrt(r_norm_2/b_norm_2)
                   << " (threshold was " << residual_tolerance_threshold << ")";
    } else {
        msg_info() << "CG diverged with a residual of |r|/|b| = " << sqrt(r_norm_2/b_norm_2)
                   << " (threshold was " << residual_tolerance_threshold << ")";
    }

    end:
    sofa::helper::AdvancedTimer::valSet("nb_iterations", static_cast<float>(iteration_number+1));
}

template <class EigenMatrix_t>
template <typename Preconditioner>
bool ConjugateGradientSolver<EigenMatrix_t>::solve(const Preconditioner & precond, const Matrix & A, const Vector & b, Vector & x) {
    // Get the method parameters
    const auto & maximum_number_of_iterations = d_maximum_number_of_iterations.getValue();
    const auto & residual_tolerance_threshold = d_residual_tolerance_threshold.getValue();
    const auto & verbose = d_verbose.getValue();

    p_squared_residuals.clear();
    p_squared_residuals.reserve(maximum_number_of_iterations);

    // Declare the method variables
    FLOATING_POINT_TYPE b_norm_2 = 0., r_norm_2 = 0.; // RHS and residual squared norms
    FLOATING_POINT_TYPE rho0 = 0., rho1 = 0.; // Temporary vectors
    FLOATING_POINT_TYPE alpha, beta; // Alpha and Beta coefficients
    FLOATING_POINT_TYPE threshold; // Residual threshold
    UNSIGNED_INTEGER_TYPE iteration_number = 0; // Current iteration number
    bool converged = false;
    UNSIGNED_INTEGER_TYPE n = A.cols();
    Vector p(n), z(n); // Search directions
    Vector r(n), q(n); // Residual
    const auto zero = (std::numeric_limits<FLOATING_POINT_TYPE>::min)(); // A numerical floating point zero

    // Make sure that the right hand side isn't zero
    b_norm_2 = b.squaredNorm();
    p_squared_initial_residual = b_norm_2;
    if (b_norm_2 < EPSILON) {
        msg_info() << "Right-hand side of the system is zero, hence x = 0.";
        x.setZero();
        goto end; // The goto is important to catch the last timer call before ending the function
    }

    // Compute the tolerance w.r.t |b| since |r|/|b| < threshold is equivalent to  r^2 < b^2 * threshold^2
    // threshold = b^2 * residual_tolerance_threshold^2
    threshold = std::max(residual_tolerance_threshold*residual_tolerance_threshold*b_norm_2, zero);

    // INITIAL RESIDUAL
    r.noalias() = b - A*x;

    // Check for initial convergence
    r_norm_2 = r.squaredNorm();
    if (r_norm_2 < threshold) {
        msg_info() << "The linear system has already reached an equilibrium state";
        msg_info() << "|r|/|b| = " << sqrt(r_norm_2/b_norm_2) << ", threshold = " << residual_tolerance_threshold;
        goto end; // The goto is important to catch the last timer call before ending the function
    }

    // Compute the initial search direction
    p = precond.solve(r);
    rho0 = r.dot(p); // |M-1 * r|^2

    // ITERATIONS
    while (not converged and iteration_number < maximum_number_of_iterations) {
        Timer::stepBegin("cg_iteration");
        // 1. Computes q(k+1) = A*p(k)
        q.noalias() = A * p;

        // 2. Computes x(k+1) and r(k+1)
        alpha = rho0 / p.dot(q); // the amount we travel on the search direction
        x += alpha * p; // Updated solution x(k+1)
        r -= alpha * q; // Updated residual r(k+1)

        // 3. Computes the new residual norm
        r_norm_2 = r.squaredNorm();
        p_squared_residuals.emplace_back(r_norm_2);

        // 4. Print information on the current iteration
        msg_info_when(verbose)  << "CG iteration #" << iteration_number+1
                                << ": |r|/|b| = "   << sqrt(r_norm_2/b_norm_2)
                                << "(threshold is " << residual_tolerance_threshold << ")";

        // 5. Check for convergence: |r|/|b| < threshold
        if (r_norm_2 < threshold) {
            converged = true;
        } else {
            // 6. Compute the next search direction
            z = precond.solve(r);  // approximately solve for "A z = r"
            rho1 = r.dot(z);
            beta = rho1 / rho0;
            p = z + beta*p;

            rho0 = rho1;
        }

        ++iteration_number;
        Timer::stepEnd("cg_iteration");
    }

    iteration_number--; // Reset to the actual index of the last iteration completed

    if (converged) {
        msg_info() << "CG converged in " << (iteration_number+1)
                   << " iterations with a residual of |r|/|b| = " << sqrt(r_norm_2/b_norm_2)
                   << " (threshold was " << residual_tolerance_threshold << ")";
    } else {
        msg_info() << "CG diverged with a residual of |r|/|b| = " << sqrt(r_norm_2/b_norm_2)
                   << " (threshold was " << residual_tolerance_threshold << ")";
    }

    end:
    sofa::helper::AdvancedTimer::valSet("nb_iterations", static_cast<float>(iteration_number+1));
    return converged;
}

template <class EigenMatrix_t>
bool ConjugateGradientSolver<EigenMatrix_t>::analyze_pattern() {
    auto A_ = this->A();
    if (not A_) {
        throw std::runtime_error("Tried to analyze an incompatible matrix (not an Eigen matrix).");
    }

    // Get the preconditioning method
    const PreconditioningMethod preconditioning_method = get_preconditioning_method_from_string(d_preconditioning_method.getValue().getSelectedItem());

    bool success = true;
    if (preconditioning_method == PreconditioningMethod::Identity || preconditioning_method == PreconditioningMethod::None) {
        p_identity.analyzePattern(A_->matrix());
        success = p_identity.info() == Eigen::Success;
    } else if (preconditioning_method == PreconditioningMethod::Diagonal) {
        p_diag.analyzePattern(A_->matrix());
        success = p_diag.info() == Eigen::Success;
#if EIGEN_VERSION_AT_LEAST(3,3,0)
    } else if (preconditioning_method == PreconditioningMethod::IncompleteCholesky) {
        p_ichol.analyzePattern(A_->matrix());
        success = p_ichol.info() == Eigen::Success;
#endif
    } else if (preconditioning_method == PreconditioningMethod::IncompleteLU) {
        p_iLU.analyzePattern(A_->matrix());
        success = p_iLU.info() == Eigen::Success;
    }

    return success;
}

template <class EigenMatrix_t>
bool ConjugateGradientSolver<EigenMatrix_t>::factorize() {
    auto A_ = this->A();
    if (not A_) {
        throw std::runtime_error("Tried to analyze an incompatible matrix (not an Eigen matrix).");
    }

    // Get the preconditioning method
    const PreconditioningMethod preconditioning_method = get_preconditioning_method_from_string(d_preconditioning_method.getValue().getSelectedItem());

    bool success = true;
    if (preconditioning_method == PreconditioningMethod::Identity || preconditioning_method == PreconditioningMethod::None) {
        p_identity.factorize(A_->matrix());
        success = p_identity.info() == Eigen::Success;
    } else if (preconditioning_method == PreconditioningMethod::Diagonal) {
        p_diag.factorize(A_->matrix());
        success = p_diag.info() == Eigen::Success;
#if EIGEN_VERSION_AT_LEAST(3,3,0)
    } else if (preconditioning_method == PreconditioningMethod::IncompleteCholesky) {
        p_ichol.factorize(A_->matrix());
        success = p_ichol.info() == Eigen::Success;
#endif
    } else if (preconditioning_method == PreconditioningMethod::IncompleteLU) {
        p_iLU.factorize(A_->matrix());
        success = p_iLU.info() == Eigen::Success;
    }

    return success;
}

template <class EigenMatrix_t>
bool ConjugateGradientSolver<EigenMatrix_t>::solve(const SofaCaribou::Algebra::BaseVector * F_, SofaCaribou::Algebra::BaseVector *X_) {
    const auto & F = dynamic_cast<const SofaCaribou::Algebra::EigenVector<Vector> *>(F_)->vector();
    auto & X = dynamic_cast<SofaCaribou::Algebra::EigenVector<Vector> *>(X_)->vector();
    const PreconditioningMethod preconditioning_method = get_preconditioning_method_from_string(d_preconditioning_method.getValue().getSelectedItem());

    bool converged = true;

    if (preconditioning_method == PreconditioningMethod::Identity || preconditioning_method == PreconditioningMethod::None) {
        converged = solve(p_identity, this->A()->matrix(), F, X);
    } else if (preconditioning_method == PreconditioningMethod::Diagonal) {
        converged = solve(p_diag, this->A()->matrix(), F, X);
#if EIGEN_VERSION_AT_LEAST(3,3,0)
    } else if (preconditioning_method == PreconditioningMethod::IncompleteCholesky) {
        converged = solve(p_ichol, this->A()->matrix(), F, X);
#endif
    } else if (preconditioning_method == PreconditioningMethod::IncompleteLU) {
        converged = solve(p_iLU, this->A()->matrix(), F, X);
    }

    return converged;
}

} // namespace SofaCaribou::solver
