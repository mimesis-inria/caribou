#pragma once

#include <Caribou/config.h>
#include <sofa/core/behavior/LinearSolver.h>
#include <SofaBaseLinearSolver/FullVector.h>
#include <SofaBaseLinearSolver/FullMatrix.h>
#include <SofaBaseLinearSolver/DefaultMultiMatrixAccessor.h>
#include <sofa/helper/OptionsGroup.h>
#include <Eigen/Core>

#include <Eigen/IterativeLinearSolvers>

namespace SofaCaribou::solver {

using sofa::core::objectmodel::Data;
using namespace sofa::component::linearsolver;
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
class ConjugateGradientSolver : public LinearSolver {

public:
    SOFA_CLASS(ConjugateGradientSolver, LinearSolver);
    using SparseMatrix = Eigen::SparseMatrix<FLOATING_POINT_TYPE, Eigen::ColMajor>;
    using Vector = Eigen::Matrix<FLOATING_POINT_TYPE, Eigen::Dynamic, 1>;

    /// Preconditioning methods
    enum class PreconditioningMethod : unsigned int {
        /// No preconditioning, hence the complete matrix won't be built. (default)
        None = 0,

        /// A naive preconditioner which approximates any matrix as the identity matrix.
        Identity = 1,

        /// Preconditioning using an approximation of A.x = b by ignoring all off-diagonal entries of A.
        Diagonal = 2,

#if EIGEN_VERSION_AT_LEAST(3,3,0)
        /// Preconditioning using an approximation of A'A.x = A'.b and ignoring all off-diagonal entries of A' A.
        LeastSquareDiagonal = 3,

        /// Preconditioning based on the incomplete Cholesky factorization.
        IncompleteCholesky = 4,
#endif

        /// Preconditioning based on the incomplete LU factorization.
        IncompleteLU = 5
    };

    /**
     * Reset the complete system (A, x and b are cleared).
     *
     * When using no preconditioner (None), this does absolutely nothing here since the complete
     * system is never built.
     *
     * This method is called by the MultiMatrix::clear() and MultiMatrix::reset() methods.
     */
    void resetSystem() final;

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
     */
    void setSystemMBKMatrix(const sofa::core::MechanicalParams* mparams) final;

    /**
     * Gives the identifier of the right-hand side vector b. This identifier will be used to find the actual vector
     * in the mechanical objects of the system.When using a preconditioner (other than None), the complete
     * dense vector is accumulated from the mechanical objects found in the graph subtree of the current context.
     */
    void setSystemRHVector(sofa::core::MultiVecDerivId b_id) final;

    /**
     * Gives the identifier of the left-hand side vector x. This identifier will be used to find the actual vector
     * in the mechanical objects of the system. When using a preconditioner (other than None), the complete
     * dense vector is accumulated from the mechanical objects found in the graph subtree of the current context.
     */
    void setSystemLHVector(sofa::core::MultiVecDerivId x_id) final;

    /**
     * Solves the system by the conjugate gradient method using the coefficients m, b and k; and
     * the vectors x and b.
     */
    void solveSystem() final;

    /**
     * List of squared residual norms (||r||^2) of every CG iterations of the last solve call.
     */
    auto squared_residuals() const -> const std::vector<FLOATING_POINT_TYPE> & {
        return p_squared_residuals;
    }

    /**
     * Squared residual norm (||r||^2) of the last right-hand side term (b in Ax=b) of the last solve call.
     */
    auto squared_initial_residual() const -> const FLOATING_POINT_TYPE & {
        return p_squared_initial_residual;
    }

    /**
     * Get the system matrix.
     */
    auto system_matrix () const -> const SparseMatrix & {
        return p_A;
    }

    /**
     * Assemble the system matrix A = (mM + bB + kK).
     * @param mparams Mechanical parameters containing the m, b and k factors.
     */
    void assemble (const sofa::core::MechanicalParams* mparams);

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
    void solve(MultiVecDeriv & b, MultiVecDeriv & x);

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
     */
    template <typename Matrix, typename Preconditioner>
    void solve(const Preconditioner & precond, const Matrix & A, const Vector & b, Vector & x);

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

    ///< Accessor used to determine the index of each mechanical object matrix and vector in the global system.
    DefaultMultiMatrixAccessor p_accessor;

    ///< The identifier of the b vector
    sofa::core::MultiVecDerivId p_b_id;

    ///< The identifier of the x vector
    sofa::core::MultiVecDerivId p_x_id;

    ///< Global system matrix (only built when a preconditioning method needs it)
    SparseMatrix p_A;

    ///< Global system solution vector (usually filled with an initial guess or the previous solution)
    Vector p_x;

    ///< Global system right-hand side vector
    Vector p_b;

    ///< Identity preconditioner
    Eigen::IdentityPreconditioner p_identity;

    ///< Diagonal preconditioner
    Eigen::DiagonalPreconditioner<FLOATING_POINT_TYPE> p_diag;

#if EIGEN_VERSION_AT_LEAST(3,3,0)
    ///< Least-Square Diagonal preconditioner
    Eigen::LeastSquareDiagonalPreconditioner<FLOATING_POINT_TYPE> p_ls_diag;

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

} // namespace SofaCaribou::solver

