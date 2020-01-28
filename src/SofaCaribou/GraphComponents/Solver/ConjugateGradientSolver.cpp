#include "ConjugateGradientSolver.h"
#include <Caribou/macros.h>
#include <sofa/core/ObjectFactory.h>
#include <sofa/helper/AdvancedTimer.h>
#include <sofa/simulation/MechanicalOperations.h>
#include <sofa/simulation/VectorOperations.h>
#include <SofaBaseLinearSolver/FullMatrix.h>
#include <SofaEigen2Solver/EigenVectorWrapper.h>
#include <iomanip>

namespace SofaCaribou::GraphComponents::solver {

using Timer = sofa::helper::AdvancedTimer;

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
, d_preconditioning_method(initData(&d_preconditioning_method,
    "preconditioning_method",
    R"(
        Preconditioning methods are:
            Identity:           No preconditioning, hence the complete matrix won't be built. (default)
            IncompleteCholesky: Preconditioning based on the incomplete Cholesky factorization.
            IncompleteLU:       Preconditioning based on the incomplete LU factorization.
    )",
    true /*displayed_in_GUI*/, false /*read_only_in_GUI*/))
{
    d_preconditioning_method.setValue(sofa::helper::OptionsGroup(std::vector<std::string> {
            "Identity",
            "IncompleteCholesky",
            "IncompleteLU"
    }));

    sofa::helper::WriteAccessor<Data< sofa::helper::OptionsGroup >> preconditioning_method = d_preconditioning_method;
    preconditioning_method->setSelectedItem((unsigned int) 0);
}

void ConjugateGradientSolver::setSystemMBKMatrix(const sofa::core::MechanicalParams* mparams) {
    Timer::stepBegin("ConjugateGradient::ComputeGlobalMatrix");
    // Save the current mechanical parameters (m, b and k factors of the mass (M), damping (B) and
    // stiffness (K) matrices)
    p_mechanical_params = mparams;

    // Get the preconditioning method
    const auto preconditioning_method = static_cast<PreconditioningMethod>(d_preconditioning_method.getValue().getSelectedId());

    // If we have a preconditioning method, the global system matrix has to be constructed from the current context subgraph.
    if (preconditioning_method != PreconditioningMethod::Identity) {
        // Step 1. Preparation stage
        //         This stage go down on the sub-graph and gather the top-level mechanical objects (mechanical objects that
        //         aren't slaved of another mechanical object through a mechanical mapping). The matrix isn't built yet,
        //         we only keep a list of mechanical objects and mappings, and get the size of global system matrix that
        //         will be built in the next stage.
        Timer::stepBegin("PrepareMatrix");

        sofa::simulation::common::MechanicalOperations mops(p_mechanical_params, this->getContext());
        p_accessor.clear();

        Timer::stepBegin("Dimension");
        mops.getMatrixDimension(nullptr, nullptr, &p_accessor);
        const auto n = static_cast<size_t>(p_accessor.getGlobalDimension());
        Timer::stepEnd("Dimension");

        Timer::stepBegin("SetupMatrixIndices");
        p_accessor.setupMatrices();
        Timer::stepEnd("SetupMatrixIndices");

        Timer::stepBegin("Clear");
        p_A.resize((int) n, (int) n);
        p_accessor.setGlobalMatrix(&p_A);
        Timer::stepEnd("Clear");

        Timer::stepEnd("PrepareMatrix");

        Timer::stepBegin("BuildMatrix");
        // Step 2. Building stage
        //         Here we go down on the current context sub-graph and call :
        //           1. ff->addMBKToMatrix(&K) for every force field "ff" found.
        //           2. pc->applyConstraint(&K) for every BaseProjectiveConstraintSet "pc" found.
        //         If a mechanical mapping "m" is found during the traversal, and m->areMatricesMapped() is false, the
        //         traversal stops in the subgraph of the mapping.
        Timer::stepBegin("TopLevelMatrices");
        mops.addMBK_ToMatrix(&p_accessor, p_mechanical_params->mFactor(), p_mechanical_params->bFactor(), p_mechanical_params->kFactor());
        Timer::stepEnd("TopLevelMatrices");
        // Step 3. Mechanical mappings
        //         In case we have mapped matrices, which is, system matrix of a slave mechanical object, accumulate its
        //         contribution to the global system matrix with:
        //           [A]ij += Jt * [A']ij * J
        //         where A is the master mechanical object's matrix, A' is the slave mechanical object matrix and J=m.getJ()
        //         is the mapping relation between the slave and its master.
        Timer::stepBegin("MappedMatrices");
        p_accessor.computeGlobalMatrix();
        Timer::stepEnd("MappedMatrices");

        // Step 4. Convert the system matrix to a compressed sparse matrix
        Timer::stepBegin("ConvertToSparse");
        p_A.compress();
        Timer::stepEnd("ConvertToSparse");
        Timer::stepEnd("BuildMatrix");

        // Step 5. Factorize the preconditioner
        Timer::stepBegin("PreconditionerFactorization");
        if (preconditioning_method == PreconditioningMethod::IncompleteCholesky) {
            p_ichol.compute(p_A.compressedMatrix);
        } else if (preconditioning_method == PreconditioningMethod::IncompleteLU) {
            p_iLU.compute(p_A.compressedMatrix);
        }
        Timer::stepEnd("PreconditionerFactorization");
    }

    Timer::stepEnd("ConjugateGradient::ComputeGlobalMatrix");
}

void ConjugateGradientSolver::setSystemRHVector(sofa::core::MultiVecDerivId b_id) {
    Timer::stepBegin("ConjugateGradient::InitialResidualVector");

    sofa::simulation::common::MechanicalOperations mop(p_mechanical_params, this->getContext());
    p_b_id = b_id;

    // Get the preconditioning method
    const auto preconditioning_method = static_cast<PreconditioningMethod>(d_preconditioning_method.getValue().getSelectedId());

    // If we have a preconditioning method, we copy the vectors of the mechanical objects into a global eigen vector.
    if (preconditioning_method != PreconditioningMethod::Identity) {
        p_b.resize(p_A.rowSize());
        EigenVectorWrapper<FLOATING_POINT_TYPE> b(p_b);
        mop.multiVector2BaseVector(p_b_id, &b, &p_accessor);
    }

    Timer::stepEnd("ConjugateGradient::InitialResidualVector");
}

void ConjugateGradientSolver::setSystemLHVector(sofa::core::MultiVecDerivId x_id) {
    Timer::stepBegin("ConjugateGradient::InitialSolutionVector");

    sofa::simulation::common::MechanicalOperations mop(p_mechanical_params, this->getContext());
    p_x_id = x_id;

    // Get the preconditioning method
    const auto preconditioning_method = static_cast<PreconditioningMethod>(d_preconditioning_method.getValue().getSelectedId());

    // If we have a preconditioning method, we copy the vectors of the mechanical objects into a global eigen vector.
    if (preconditioning_method != PreconditioningMethod::Identity) {
        p_x.resize(p_A.rowSize());
        EigenVectorWrapper<FLOATING_POINT_TYPE> x(p_x);
        mop.multiVector2BaseVector(p_x_id, &x, &p_accessor);
    }

    Timer::stepEnd("ConjugateGradient::InitialSolutionVector");
}

void ConjugateGradientSolver::solve(MultiVecDeriv & b, MultiVecDeriv & x) {
    sofa::simulation::common::VectorOperations vop( p_mechanical_params, this->getContext() );
    sofa::simulation::common::MechanicalOperations mop( p_mechanical_params, this->getContext() );

    // Create temporary vectors needed for the method
    MultiVecDeriv r(&vop);
    MultiVecDeriv p(&vop);
    MultiVecDeriv q(&vop);

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

    // Make sure that the right hand side isn't zero
    b_norm = b.norm();
    if (IN_OPEN_INTERVAL(-EPSILON, b_norm, EPSILON)) {
        msg_info() << "Right-hand side of the system is zero, hence x = 0.";
        x.clear();
        return;
    }

    // INITIAL RESIDUAL
    // Do the A*x(0) with visitors since we did not construct the matrix A
    mop.propagateDxAndResetDf(x, q); // Set q = 0 and calls applyJ(x) on every mechanical mappings
    mop.addMBKdx(q, m_coef, b_coef, k_coef, false); // q = (m M + b B + k K) x

    // Do the projection of the result in the constrained space
    mop.projectResponse(q); // BaseProjectiveConstraintSet::projectResponse(q)

    // Finally, compute the initial residual r = b - A*x
    r.eq( b, q, -1.0 );   // r = b - q

    // Check for initial convergence: |r0|/|b| < threshold
    rho0 = r.dot(r);
    r_norm = sqrt(rho0);
    if (r_norm < residual_tolerance_threshold*b_norm) {
        msg_info() << "The linear system has already reached an equilibrium state";
        msg_info() << "|R| = " << r_norm << ", |b| = " << b_norm << ", threshold = " << residual_tolerance_threshold;
        return;
    }

    // ITERATIONS
    p = r; // p(0) = r(0)
    while (iteration_number < maximum_number_of_iterations) {
        // 1. Computes q(k+1) = A*p(k)
        mop.propagateDxAndResetDf(p, q); // Set q = 0 and calls applyJ(p) on every mechanical mappings
        mop.addMBKdx(q, m_coef, b_coef, k_coef, false); // q = (m M + b B + k K) x
        mop.projectResponse(q); // BaseProjectiveConstraintSet::projectResponse(q)

        // 2. Computes x(k+1) and r(k+1)
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

template <typename Matrix, typename Preconditioner>
void ConjugateGradientSolver::solve(const Preconditioner & precond, const Matrix & A, const Vector & b, Vector & x) {
    // Get the method parameters
    const auto & maximum_number_of_iterations = d_maximum_number_of_iterations.getValue();
    const auto & residual_tolerance_threshold = d_residual_tolerance_threshold.getValue();
    const auto zero = (std::numeric_limits<FLOATING_POINT_TYPE>::min)();

    // Declare the method variables
    FLOATING_POINT_TYPE b_norm_2; // RHS norm
    FLOATING_POINT_TYPE rho0, rho1; // Stores r*r as it is used two times per iterations
    FLOATING_POINT_TYPE alpha, beta; // Alpha and Beta coefficients
    UNSIGNED_INTEGER_TYPE iteration_number = 0; // Current iteration number
    UNSIGNED_INTEGER_TYPE n = A.cols();

    // Make sure that the right hand side isn't zero
    b_norm_2 = b.squaredNorm();
    if (b_norm_2 < EPSILON) {
        msg_info() << "Right-hand side of the system is zero, hence x = 0.";
        x.setZero();
        return;
    }

    // Compute the tolerance w.r.t |b| since |r|/|b| < threshold is equivalent to  r^2 < b^2 * threshold^2
    // threshold = b^2 * residual_tolerance_threshold^2
    auto threshold = Eigen::numext::maxi(residual_tolerance_threshold*residual_tolerance_threshold*b_norm_2,zero);

    // INITIAL RESIDUAL
    Vector r = b - A*x;

    // Check for initial convergence
    rho0 = r.dot(r);
    if (rho0 < threshold) {
        msg_info() << "The linear system has already reached an equilibrium state";
        msg_info() << "|r0|/|b| = " << sqrt(rho0/b_norm_2) << ", threshold = " << residual_tolerance_threshold;
        return;
    }

    // Compute the initial search direction
    Vector p(n), z(n);
    p = precond.solve(r);
    rho0 = r.dot(p); // |M-1 * r|^2

    // ITERATIONS
    while (iteration_number < maximum_number_of_iterations) {
        // 1. Computes q(k+1) = A*p(k)
        Vector q = A * p;

        // 2. Computes x(k+1) and r(k+1)
        alpha = rho0 / p.dot(q); // the amount we travel on the search direction
        x += alpha * p; // Updated solution x(k+1)
        r -= alpha * q; // Updated residual r(k+1)

        // 3. Computes the new residual norm
        rho1 = r.dot(r);

        // 4. Print information on the current iteration
        msg_info()  << "CG iteration #" << iteration_number+1
                    << ": |r0|/|b| = "  << sqrt(rho1/b_norm_2)
                    << ", threshold = " << residual_tolerance_threshold;

        // 5. Check for convergence: |r|/|b| < threshold
        if (rho0 < threshold) {
            msg_info() << "CG converged!";
            return;
        }

        // 6. Compute the next search direction
        z = precond.solve(r);  // approximately solve for "A z = r"
        rho1 = r.dot(z);
        beta = rho1 / rho0;
        p = z + beta*p;

        rho0 = rho1;
        ++iteration_number;
    }
}

void ConjugateGradientSolver::solveSystem() {
    sofa::simulation::common::VectorOperations vop( p_mechanical_params, this->getContext() );
    sofa::simulation::common::MechanicalOperations mop( p_mechanical_params, this->getContext() );

    // Gather the x and b vector identifiers
    MultiVecDeriv x(&vop, p_x_id);
    MultiVecDeriv b(&vop, p_b_id);

    // Get the preconditioning method
    const auto preconditioning_method = static_cast<PreconditioningMethod>(d_preconditioning_method.getValue().getSelectedId());


    Timer::stepBegin("ConjugateGradient::solve");
    if (preconditioning_method == PreconditioningMethod::Identity) {
        // Solve without having filled the global matrix A (not needed since no preconditioning)
        solve(b, x);
    } else {
        // Solve using a preconditioning method. Here the global matrix A and the vectors x and b have been built
        // previously during the calls to setSystemMBKMatrix, setSystemLHVector and setSystemRHVector, respectively.

        if (preconditioning_method == PreconditioningMethod::IncompleteCholesky) {
            solve(p_ichol, p_A.compressedMatrix, p_b, p_x);
        } else if (preconditioning_method == PreconditioningMethod::IncompleteLU) {
            solve(p_iLU, p_A.compressedMatrix, p_b, p_x);
        }

        // Copy the solution into the mechanical objects of the current context sub-graph.
        EigenVectorWrapper<FLOATING_POINT_TYPE> x_wrapper(p_x);
        mop.baseVector2MultiVector(p_x_id, &x_wrapper, &p_accessor);
    }

    Timer::stepEnd("ConjugateGradient::solve");
}

int CGLinearSolverClass = sofa::core::RegisterObject("Linear system solver using the conjugate gradient iterative algorithm")
                              .add< ConjugateGradientSolver>(true);

} // namespace SofaCaribou::GraphComponents::solver