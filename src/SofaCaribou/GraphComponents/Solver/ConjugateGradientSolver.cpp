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
    )",
    true /*displayed_in_GUI*/, false /*read_only_in_GUI*/))
{
    d_preconditioning_method.setValue(sofa::helper::OptionsGroup(std::vector<std::string> {
            "Identity",
            "IncompleteCholesky"
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

    // If no preconditioning, we are done.
    if (preconditioning_method == PreconditioningMethod::Identity) {
        return;
    }

    // Else, the global system matrix has to be constructed from the current context subgraph.

    // Step 1. Preparation stage
    //         This stage go down on the sub-graph and gather the top-level mechanical objects (mechanical objects that
    //         aren't slaved of another mechanical object through a mechanical mapping). The matrix isn't built yet,
    //         we only keep a list of mechanical objects and mappings, and get the size of global system matrix that
    //         will be built in the next stage.
    Timer::stepBegin("PrepareMatrix");

    sofa::simulation::common::MechanicalOperations mops(p_mechanical_params, this->getContext());
    FullMatrix<FLOATING_POINT_TYPE> A;
    p_accessor.setGlobalMatrix(&A);
    p_accessor.clear();

    Timer::stepBegin("Dimension");
    mops.getMatrixDimension(nullptr, nullptr, &p_accessor);
    const auto n = static_cast<size_t>(p_accessor.getGlobalDimension());
    Timer::stepEnd("Dimension");

    Timer::stepBegin("SetupMatrixIndices");
    p_accessor.setupMatrices();
    Timer::stepEnd("SetupMatrixIndices");

    Timer::stepBegin("Clear");
    A.resize((int) n, (int) n); // Resize and sets all entries to 0
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
    std::vector<Eigen::Triplet<double>> entries;
    entries.reserve(n*3);
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            const auto & v = A[i][j];
            if (not IN_CLOSED_INTERVAL(-EPSILON, v, EPSILON)) {
                entries.emplace_back(i, j, v);
            }
        }
    }
    entries.shrink_to_fit();
    p_A.resize(n, n);
    p_A.setFromTriplets(entries.begin(), entries.end());
    p_A.makeCompressed();
    Timer::stepEnd("ConvertToSparse");
    Timer::stepEnd("BuildMatrix");

    Timer::stepEnd("ConjugateGradient::ComputeGlobalMatrix");
}

void ConjugateGradientSolver::setSystemRHVector(sofa::core::MultiVecDerivId b_id) {
    Timer::stepBegin("ConjugateGradient::InitialResidualVector");

    sofa::simulation::common::MechanicalOperations mop(p_mechanical_params, this->getContext());
    p_b_id = b_id;

    // Get the preconditioning method
    const auto preconditioning_method = static_cast<PreconditioningMethod>(d_preconditioning_method.getValue().getSelectedId());

    // If no preconditioning, we are done.
    if (preconditioning_method == PreconditioningMethod::Identity) {
        return;
    }

    p_b.resize(p_A.rows());
    EigenVectorWrapper<FLOATING_POINT_TYPE> b(p_b);
    mop.multiVector2BaseVector(p_b_id, &b, &p_accessor);

    Timer::stepEnd("ConjugateGradient::InitialResidualVector");
}

void ConjugateGradientSolver::setSystemLHVector(sofa::core::MultiVecDerivId x_id) {
    Timer::stepBegin("ConjugateGradient::InitialSolutionVector");

    sofa::simulation::common::MechanicalOperations mop(p_mechanical_params, this->getContext());
    p_x_id = x_id;

    // Get the preconditioning method
    const auto preconditioning_method = static_cast<PreconditioningMethod>(d_preconditioning_method.getValue().getSelectedId());

    // If no preconditioning, we are done.
    if (preconditioning_method == PreconditioningMethod::Identity) {
        return;
    }

    p_x.resize(p_A.rows());
    EigenVectorWrapper<FLOATING_POINT_TYPE> x(p_x);
    mop.multiVector2BaseVector(p_x_id, &x, &p_accessor);

    Timer::stepEnd("ConjugateGradient::InitialSolutionVector");
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

    Timer::stepBegin("ConjugateGradient::solve");

    // Compute the first residual r(0) = b - A*x(0) = b - q
    mop.propagateDxAndResetDf(x, q); // Set q = 0 and calls applyJ(x) on every mechanical mappings
    mop.addMBKdx(q, m_coef, b_coef, k_coef, false); // q = (m M + b B + k K) x
    mop.projectResponse(q); // BaseProjectiveConstraintSet::projectResponse(q)
    r.eq( b, q, -1.0 );   // r = b - q

    {
        Eigen::Matrix<FLOATING_POINT_TYPE, Eigen::Dynamic, 1> rr = (p_b - (p_A * p_x));
        const auto n = p_A.rows();
        std::cout << "[ ";
        for (size_t i = 0; i < n; ++i)
            std::cout << std::setprecision(3) << rr[i] << " ";
        EigenVectorWrapper<FLOATING_POINT_TYPE> rrr(rr);
        std::cout<<std::endl<<r<<std::endl;
        mop.multiVectorPeqBaseVector(p.id(), &rrr, &p_accessor);
        mop.projectResponse(p);
//        std::cout<<std::endl<<p<<std::endl;
        p.peq(r, -1.0);
        std::cout<< "!_!_!_!_!_!_! DIFF IS " << sqrt(p.dot(p)) << std::endl;
    }

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

    Timer::stepEnd("ConjugateGradient::solve");
}

int CGLinearSolverClass = sofa::core::RegisterObject("Linear system solver using the conjugate gradient iterative algorithm")
                              .add< ConjugateGradientSolver>(true);

} // namespace SofaCaribou::GraphComponents::solver