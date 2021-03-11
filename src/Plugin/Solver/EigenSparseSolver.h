#pragma once

#include <Caribou/config.h>
#include <SofaCaribou/Algebra/EigenMatrix.h>
#include <SofaCaribou/Algebra/EigenVector.h>
#include <SofaCaribou/Solver/LinearSolver.h>

#include <sofa/core/behavior/LinearSolver.h>
#include <SofaBaseLinearSolver/DefaultMultiMatrixAccessor.h>
#include <Eigen/Sparse>

namespace SofaCaribou::solver {

namespace internal {
template <typename T>
struct solver_traits {};
}

/**
 * Base class for sparse direct solvers using Eigen as a backend solver.
 *
 * Note that the complete system matrix has to be assemble for any of the derivated solvers.
 * This means that the mass, damping and forcefield components in the scene graph need to implement
 * the addMtoMatrix, addBtoMatrix and addMtoMatrix methods respectively.
 *
 * @tparam EigenSolver_t Eigen solver type (eg.: SimplicialLLT, SparseLU, PardisoLLT, etc.)
 */
template <class EigenSolver_t>
class EigenSparseSolver : public sofa::core::behavior::LinearSolver, public SofaCaribou::solver::LinearSolver {
public:
    SOFA_CLASS(SOFA_TEMPLATE(EigenSparseSolver, EigenSolver_t), sofa::core::behavior::LinearSolver);
    using EigenSolver = EigenSolver_t;
    using SparseMatrix = std::remove_cv_t<typename EigenSolver::MatrixType>;
    using Vector = Eigen::Matrix<FLOATING_POINT_TYPE, Eigen::Dynamic, 1>;

    /**
     * Assemble the system matrix A = (mM + bB + kK) inside the SparseMatrix p_A.
     * @param mparams Mechanical parameters containing the m, b and k factors.
     */
    virtual void assemble (const sofa::core::MechanicalParams* mparams);

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
     * @param mparams Contains the coefficients m, b and k of the matrices M, B and K
     */
    void setSystemMBKMatrix(const sofa::core::MechanicalParams* mparams) final;

    /**
     * Gives the identifier of the right-hand side vector b. This identifier will be used to find the actual vector
     * in the mechanical objects of the system. The complete dense vector is accumulated from the mechanical objects
     * found in the graph subtree of the current context.
     */
    void setSystemRHVector(sofa::core::MultiVecDerivId b_id) final;

    /**
     * Gives the identifier of the left-hand side vector x. This identifier will be used to find the actual vector
     * in the mechanical objects of the system. The complete dense vector is accumulated from the mechanical objects
     * found in the graph subtree of the current context.
     */
    void setSystemLHVector(sofa::core::MultiVecDerivId x_id) final;

    /** Solves the system using the Eigen solver. */
    void solveSystem() override;

    /**
     * States if the system matrix is symmetric. Note that this value isn't set automatically, the user must
     * explicitly specify it using set_symmetric(true). When it is true, some optimizations will be enabled.
     */
    inline virtual auto symmetric() const -> bool {return p_is_symmetric;}

    /** Explicitly states if this matrix is symmetric. */
    inline virtual void set_symmetric(bool is_symmetric) { p_is_symmetric = is_symmetric; }

    /** Get a readonly reference to the backend solver */
    auto solver() const -> const EigenSolver & { return p_solver; }

    /** Get a reference to the backend solver */
    auto solver() -> EigenSolver & { return p_solver; }

    /** Get a readonly reference to the mechanical parameters */
    auto mechanical_params() const -> const sofa::core::MechanicalParams & { return p_mechanical_params; }


    /**
     * Get a readonly reference to the multi-matrix accessor. This accessor doesn't actually hold the system matrix, but
     * will be the one responsible for accumulating it. It can also be use to compute the stiffness or the mass matrix
     * of a given mechanical node or mapped node.
     *
     * Use the EigenSparseSolver::matrix() method to get a reference to the global system matrix.
     */
    auto matrix_accessor() const -> const sofa::component::linearsolver::DefaultMultiMatrixAccessor & { return p_accessor; }

    /** Get a readonly reference to the global assembled system matrix. */
    auto A() const -> const SparseMatrix & { return p_A; }

    /** Get a readonly reference to the right-hand side vector identifier */
    auto b_id() const -> const sofa::core::MultiVecDerivId & { return p_b_id; }

    /** Get a readonly reference to the left-hand side unknown vector identifier */
    auto x_id() const -> const sofa::core::MultiVecDerivId & { return p_x_id; }

    /** True if the solver has successfully factorize the system matrix */
    auto A_is_factorized() const -> bool { return p_A_is_factorized; }

    // SOFA overrides
    static auto GetCustomTemplateName() -> std::string;

    template<class Derived>
    static auto canCreate(Derived* o, sofa::core::objectmodel::BaseContext* context, sofa::core::objectmodel::BaseObjectDescription* arg) -> bool;

private:
    /**
     * @see SofaCaribou::solver::LinearSolver::create_new_matrix
     */
    sofa::defaulttype::BaseMatrix * create_new_matrix(sofa::Size rows, sofa::Size cols) const override {
        auto * matrix = new SofaCaribou::Algebra::EigenMatrix<SparseMatrix> (static_cast<Eigen::Index>(rows), static_cast<Eigen::Index>(cols));
        if (symmetric()) {
            matrix->set_symmetric(symmetric());
        }
        return matrix;
    }

    /**
     * @see SofaCaribou::solver::LinearSolver::create_new_vector
     */
    sofa::defaulttype::BaseVector * create_new_vector(sofa::Size n) const override {
        return new SofaCaribou::Algebra::EigenVector<Vector>(n);
    }

    /**
     * @see SofaCaribou::solver::LinearSolver::solve
     */
    bool solve(const sofa::defaulttype::BaseMatrix * A,
               const sofa::defaulttype::BaseVector * F,
               sofa::defaulttype::BaseVector * X) const override;

private:
    /// Private members

    ///< The Eigen solver used (its type is passed as a template parameter and must be derived from Eigen::SparseSolverBase)
    EigenSolver p_solver;

    ///< The mechanical parameters containing the m, b and k coefficients.
    sofa::core::MechanicalParams p_mechanical_params;

    ///< Accessor used to determine the index of each mechanical object matrix and vector in the global system.
    sofa::component::linearsolver::DefaultMultiMatrixAccessor p_accessor;

    ///< The identifier of the b vector
    sofa::core::MultiVecDerivId p_b_id;

    ///< The identifier of the x vector
    sofa::core::MultiVecDerivId p_x_id;

    ///< Global system matrix
    SparseMatrix p_A;

    ///< Global system solution vector (usually filled with an initial guess or the previous solution)
    Vector p_x;

    ///< Global system right-hand side vector
    Vector p_b;

    ///< True if the solver has successfully factorize the system matrix
    bool p_A_is_factorized;

    ///< States if the system matrix is symmetric. Note that this value isn't set automatically, the user must
    ///< explicitly specify it using set_symmetric(true). When it is true, some optimizations will be enabled.
    bool p_is_symmetric = false;
};

} // namespace SofaCaribou::solver
