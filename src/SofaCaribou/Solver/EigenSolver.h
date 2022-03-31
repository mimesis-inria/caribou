#pragma once

#include <SofaCaribou/config.h>
#include <SofaCaribou/Algebra/EigenMatrix.h>
#include <SofaCaribou/Algebra/EigenVector.h>
#include <SofaCaribou/Solver/LinearSolver.h>

DISABLE_ALL_WARNINGS_BEGIN
#include <sofa/version.h>
#include <sofa/core/MechanicalParams.h>
#include <sofa/core/behavior/LinearSolver.h>
#include <SofaBaseLinearSolver/DefaultMultiMatrixAccessor.h>
DISABLE_ALL_WARNINGS_END

#if (defined(SOFA_VERSION) && SOFA_VERSION < 201200)
namespace sofa {
using Size = unsigned int;
using Index = unsigned int;
}
#endif

namespace SofaCaribou::solver {

/**
 * Base class for linear solvers using Eigen matrices.
 *
 * Note that the complete system matrix has to be assemble for any of the derived solvers.
 * This means that the mass, damping and forcefield components in the scene graph need to implement
 * the addMtoMatrix, addBtoMatrix and addMtoMatrix methods respectively.
 *
 * @tparam EigenMatrix_t Eigen matrix type (eg.: SparseMatrix, Matrix, DiagonalMatrix, etc.)
 */
template <class EigenMatrix_t>
class EigenSolver : public sofa::core::behavior::LinearSolver, public SofaCaribou::solver::LinearSolver {
public:
    SOFA_CLASS(SOFA_TEMPLATE(EigenSolver, EigenMatrix_t), sofa::core::behavior::LinearSolver);
    using Scalar = typename Eigen::MatrixBase<EigenMatrix_t>::Scalar;
    using Matrix = EigenMatrix_t;
    using Vector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;

    EigenSolver() = default;

    /**
     * Assemble the system matrix A = (mM + bB + kK).
     *
     * @param mparams Mechanical parameters containing the m, b and k factors.
     * @return The matrix accessor containing the lists of top level mechanical objects, a pointer for
     * their matrix and a vector of mappings for mapped mechanical objects.
     */
    auto assemble (const sofa::core::MechanicalParams* mparams) -> sofa::component::linearsolver::DefaultMultiMatrixAccessor;

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
    void setSystemMBKMatrix(const sofa::core::MechanicalParams* mparams) override;

    /**
     * Gives the identifier of the right-hand side vector b. This identifier will be used to find the actual vector
     * in the mechanical objects of the system. The complete dense vector is accumulated from the mechanical objects
     * found in the graph subtree of the current context.
     */
    void setSystemRHVector(sofa::core::MultiVecDerivId b_id) override;

    /**
     * Gives the identifier of the left-hand side vector x. This identifier will be used to find the actual vector
     * in the mechanical objects of the system. The complete dense vector is accumulated from the mechanical objects
     * found in the graph subtree of the current context.
     */
    void setSystemLHVector(sofa::core::MultiVecDerivId x_id) override;

    /** Solves the system using the Eigen solver. */
    void solveSystem() override;

    /**
     * States if the system matrix is symmetric. Note that this value isn't set automatically, the user must
     * explicitly specify it using set_symmetric(true). When it is true, some optimizations will be enabled.
     */
    inline virtual auto symmetric() const -> bool {return p_is_symmetric;}

    /** Explicitly states if this matrix is symmetric. */
    inline virtual void set_symmetric(bool is_symmetric) { p_is_symmetric = is_symmetric; }

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
    auto A() const -> const SofaCaribou::Algebra::EigenMatrix<Matrix> * { return p_A_ptr; }

    /** Get a readonly reference to the right-hand side vector. */
    auto b() const -> const SofaCaribou::Algebra::EigenVector<Vector> * { return &p_b;}

    /** Get a readonly reference to the right-hand side vector identifier */
    auto b_id() const -> const sofa::core::MultiVecDerivId & { return p_b_id; }

    /** Get a readonly reference to the left-hand side unknown vector. */
    auto x() const -> const SofaCaribou::Algebra::EigenVector<Vector> * { return &p_x;}

    /** Get a readonly reference to the left-hand side unknown vector identifier */
    auto x_id() const -> const sofa::core::MultiVecDerivId & { return p_x_id; }

    /** True if the solver has successfully factorize the system matrix */
    auto A_is_factorized() const -> bool { return p_A_is_factorized; }

    /** @see SofaCaribou::solver::LinearSolver::is_iterative */
    bool is_iterative() const override {
        return false;
    }

    /** @see SofaCaribou::solver::LinearSolver::squared_residuals */
    auto squared_residuals() const -> std::vector<FLOATING_POINT_TYPE> override {
        // Return an empty sets by default as most solvers inheriting from this class
        // are direct solvers.
        return {};
    }

    // SOFA overrides
    static auto GetCustomTemplateName() -> std::string;

    template<class Derived>
    static auto canCreate(Derived* o, sofa::core::objectmodel::BaseContext* context, sofa::core::objectmodel::BaseObjectDescription* arg) -> bool;

protected:
    void set_system_matrix(const sofa::defaulttype::BaseMatrix * A) override {
        p_A_ptr = dynamic_cast<const SofaCaribou::Algebra::EigenMatrix<Matrix> *>(A);
    }
private:
    /**
     * @see SofaCaribou::solver::LinearSolver::create_new_matrix
     */
    sofa::defaulttype::BaseMatrix * create_new_matrix(sofa::Size rows, sofa::Size cols) const override {
        auto * matrix = new SofaCaribou::Algebra::EigenMatrix<Matrix> (static_cast<Eigen::Index>(rows), static_cast<Eigen::Index>(cols));
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

private:
    /// Private members

    /// The mechanical parameters containing the m, b and k coefficients.
    sofa::core::MechanicalParams p_mechanical_params;

    /// Accessor used to determine the index of each mechanical object matrix and vector in the global system.
    sofa::component::linearsolver::DefaultMultiMatrixAccessor p_accessor;

    /// The identifier of the b vector
    sofa::core::MultiVecDerivId p_b_id;

    /// The identifier of the x vector
    sofa::core::MultiVecDerivId p_x_id;

    /// Global system matrix (assembled using legacy ODE solvers)
    SofaCaribou::Algebra::EigenMatrix<Matrix> p_A;

    /// Pointer to the global system matrix. For legacy ODE solvers, it points to
    /// p_A which is assembled by this linear solver. For Caribou's ODE solvers, it
    /// points to the system matrix assembled by the ODE itself.
    const SofaCaribou::Algebra::EigenMatrix<Matrix> * p_A_ptr;

    /// Global system solution vector (usually filled with an initial guess or the previous solution)
    SofaCaribou::Algebra::EigenVector<Vector> p_x;

    /// Global system right-hand side vector
    SofaCaribou::Algebra::EigenVector<Vector> p_b;

    /// True if the solver has successfully factorize the system matrix
    bool p_A_is_factorized {};

    /// States if the system matrix is symmetric. Note that this value isn't set automatically, the user must
    /// explicitly specify it using set_symmetric(true). When it is true, some optimizations will be enabled.
    bool p_is_symmetric = false;
};

extern template class EigenSolver<Eigen::SparseMatrix<FLOATING_POINT_TYPE, Eigen::ColMajor, int>>;
extern template class EigenSolver<Eigen::SparseMatrix<FLOATING_POINT_TYPE, Eigen::RowMajor, int>>;

} // namespace SofaCaribou::solver
