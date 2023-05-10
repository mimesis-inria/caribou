#pragma once

#include <SofaCaribou/config.h>
#include <SofaCaribou/Solver/EigenSolver.h>

DISABLE_ALL_WARNINGS_BEGIN
#include <sofa/core/behavior/MultiVec.h>
#include <sofa/core/MechanicalParams.h>
#include <sofa/helper/OptionsGroup.h>
DISABLE_ALL_WARNINGS_END

#include <Eigen/Core>
#include <Eigen/IterativeLinearSolvers>

namespace SofaCaribou::solver {

using sofa::core::objectmodel::Data;
using namespace sofa::core::behavior;

/**
 * Implements a Conjugate Gradient method for solving a linear system of the type Ax = b.
 *
 * When using no preconditioner, no actual system is built by this component, meaning that
 * the matrix A, and the vectors x and b are never accumulated into memory. Instead, this
 * solver relies on the vectors stored in the mechanical objects of the current scene context
 * graph. The multiplication of the matrix A and a given vector x is computed by accumulating
 * the result b=Ax from ff.addDForce(x, b) of every forcefields acting on a given mechanical
 * object.
 *
 * When using a preconditioner, the matrix has to be assembled since the preconditioner needs
 * to factorize it. In this case, the complete system matrix A and dense vector b are first
 * accumulated from the mechanical objects of the current scene context graph. Once the dense
 * vector x is found, it is propagated back to the mechanical object's vectors.
 */
template <class EigenMatrix_t>
class ConjugateGradientSolver : public EigenSolver<EigenMatrix_t> {

public:
    SOFA_CLASS(SOFA_TEMPLATE(ConjugateGradientSolver, EigenMatrix_t), SOFA_TEMPLATE(EigenSolver, EigenMatrix_t));
    using Base = EigenSolver<EigenMatrix_t>;
    using Matrix = typename Base::Matrix;
    using Vector = typename Base::Vector;

    /// Preconditioning methods
    enum class PreconditioningMethod : unsigned int {
        /// No preconditioning, hence the complete matrix won't be built. (default)
        None = 0,

        /// A naive preconditioner which approximates any matrix as the identity matrix.
        Identity = 1,

        /// Preconditioning using an approximation of A.x = b by ignoring all off-diagonal entries of A.
        Diagonal = 2,

#if EIGEN_VERSION_AT_LEAST(3,3,0)
        /// Preconditioning based on the incomplete Cholesky factorization.
        IncompleteCholesky = 4,
#endif

        /// Preconditioning based on the incomplete LU factorization.
        IncompleteLU = 5
    };

    /**
     * Set the linear system matrix A = (mM + bB + kK), storing the coefficients m, b and k of
     * the mechanical M,B,K matrices.
     *
     * When using no preconditioner (None), no actual matrix is built by this method since it is not
     * required by the Conjugate Gradient algorithm. Only the coefficients m, b and k are stored.
     *
     * When using a preconditioner, the complete system A is accumulated into a sparse matrix, and the
     * preconditioner factorize this resulting matrix for later use during the solve step.
     *
     * @param mparams Contains the coefficients m, b and k of the matrices M, B and K
     *
     * \note Only used by SOFA ODE solvers and the Caribou's LegacyStaticODESolver.
     */
    
    void setSystemMBKMatrix(const sofa::core::MechanicalParams* mparams) final;

    /**
     * Gives the identifier of the right-hand side vector b. This identifier will be used to find the actual vector
     * in the mechanical objects of the system.When using a preconditioner (other than None), the complete
     * dense vector is accumulated from the mechanical objects found in the graph subtree of the current context.
     *
     * \note Only used by SOFA ODE solvers and the Caribou's LegacyStaticODESolver.
     */
    
    void setSystemRHVector(sofa::core::MultiVecDerivId b_id) final;

    /**
     * Gives the identifier of the left-hand side vector x. This identifier will be used to find the actual vector
     * in the mechanical objects of the system. When using a preconditioner (other than None), the complete
     * dense vector is accumulated from the mechanical objects found in the graph subtree of the current context.
     *
     * \note Only used by SOFA ODE solvers and the Caribou's LegacyStaticODESolver.
     */
    
    void setSystemLHVector(sofa::core::MultiVecDerivId x_id) final;

    /**
     * Solves the system by the conjugate gradient method using the coefficients m, b and k; and
     * the vectors x and b.
     *
     * \note Only used by SOFA ODE solvers and the Caribou's LegacyStaticODESolver.
     */
    
    void solveSystem() final;

    /** @see SofaCaribou::solver::LinearSolver::is_iterative */
    
    bool is_iterative() const override {
        return true;
    }

    /**
     * List of squared residual norms (||r||^2) of every CG iterations of the last solve call.
     */
    std::vector<FLOATING_POINT_TYPE> squared_residuals() const override {
        return p_squared_residuals;
    }

    /**
     * Squared residual norm (||r||^2) of the last right-hand side term (b in Ax=b) of the last solve call.
     */
    auto squared_initial_residual() const -> const FLOATING_POINT_TYPE & {
        return p_squared_initial_residual;
    }

    template<typename Derived>
    static auto canCreate(Derived*, sofa::core::objectmodel::BaseContext*, sofa::core::objectmodel::BaseObjectDescription*) -> bool {
        return true;
    }
protected:
    /// Constructor
    ConjugateGradientSolver();

    /// Methods
    /**
     * Solve the linear system Ax = b using Sofa's graph scene.
     *
     * Here, the matrix A is never built. The matrix-vector product Ax is done by going down in the subgraph
     * of the current context and doing the multiplication directly in the vectors of mechanical objects found.
     *
     * \warning Since the matrix A is never built, no preconditioning method can be used with this method.
     *
     * @param b The right-hand side vector of the system
     * @param x The solution vector of the system. It should be filled with an initial guess or the previous solution.
     */
    void solve(sofa::core::behavior::MultiVecDeriv & b, sofa::core::behavior::MultiVecDeriv & x);

    /** @see LinearSolver::analyze_pattern */
    
    bool analyze_pattern() override;

    /** @see LinearSolver::factorize */
    
    bool factorize() override;

    /** @see SofaCaribou::solver::LinearSolver::solve */
    
    bool solve(const SofaCaribou::Algebra::BaseVector * F, SofaCaribou::Algebra::BaseVector * X) override;

    /**
     * Solve the linear system Ax = b using a preconditioner.
     *
     * @tparam Derived The type of the matrix A, can be a dense or a sparse matrix.
     * @tparam Preconditioner The type of the preconditioner.
     *
     * @param precond The preconditioner (must be a
     * @param A The system matrix as an Eigen matrix
     * @param b The right-hand side vector of the system
     * @param x The solution vector of the system. It should be filled with an initial guess or the previous solution.
     * @return True if the CG converged, false otherwise.
     */
    template <typename Preconditioner>
    bool solve(const Preconditioner & precond, const Matrix & A, const Vector & b, Vector & x);

    /// INPUTS
    Data<bool> d_verbose;
    Data<unsigned int> d_maximum_number_of_iterations;
    Data<FLOATING_POINT_TYPE> d_residual_tolerance_threshold;
    Data< sofa::helper::OptionsGroup > d_preconditioning_method;

private:
    /// Private methods
    /**
     * @brief Get the identifier of the preconditioning method from its name
     */
    PreconditioningMethod get_preconditioning_method_from_string(const std::string & preconditioner_name) const;

    /// Private members
    ///< The mechanical parameters containing the m, b and k coefficients.
    sofa::core::MechanicalParams p_mechanical_params;

    ///< The identifier of the b vector
    sofa::core::MultiVecDerivId p_b_id;

    ///< The identifier of the x vector
    sofa::core::MultiVecDerivId p_x_id;

    ///< Identity preconditioner
    Eigen::IdentityPreconditioner p_identity;

    ///< Diagonal preconditioner
    Eigen::DiagonalPreconditioner<FLOATING_POINT_TYPE> p_diag;

#if EIGEN_VERSION_AT_LEAST(3,3,0)
    ///< Incomplete Cholesky preconditioner
    Eigen::IncompleteCholesky<FLOATING_POINT_TYPE> p_ichol;
#endif

    ///< Incomplete LU preconditioner
    Eigen::IncompleteLUT<FLOATING_POINT_TYPE> p_iLU;

    ///< Contains the list of available preconditioners with their respective identifier
    std::vector<std::pair<std::string, PreconditioningMethod>> p_preconditioners;

    ///< List of squared residual norms (||r||^2) of every CG iterations of the last solve call.
    std::vector<FLOATING_POINT_TYPE> p_squared_residuals;

    ///< Squared residual norm (||r||^2) of the last right-hand side term (b in Ax=b) of the last solve call.
    FLOATING_POINT_TYPE p_squared_initial_residual;
};

extern template class ConjugateGradientSolver<Eigen::SparseMatrix<FLOATING_POINT_TYPE, Eigen::RowMajor, int>>;
} // namespace SofaCaribou::solver

