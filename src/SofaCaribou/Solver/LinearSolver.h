#pragma once

#include <SofaCaribou/config.h>

#if (defined(SOFA_VERSION) && SOFA_VERSION < 211200)
namespace sofa::defaulttype {
    class BaseMatrix;
    class BaseVector;
}
#else
namespace sofa::linearalgebra {
    class BaseVector;
    class BaseMatrix;
}
#endif // (defined(SOFA_VERSION) && SOFA_VERSION < 211299)

namespace SofaCaribou::solver {

/**
 * Base interface for linear solvers. This interface define the generic API that linear solvers in Caribou must
 * provide.
 */
class LinearSolver {
public:
#if (defined(SOFA_VERSION) && SOFA_VERSION < 211200)
    using BaseMatrix = sofa::defaulttype::BaseMatrix;
    using BaseVector = sofa::defaulttype::BaseVector;
#else
    using BaseMatrix = sofa::linearalgebra::BaseMatrix;
    using BaseVector = sofa::linearalgebra::BaseVector;
#endif
    /**
     * Creates a new BaseMatrix of size rows x cols.
     *
     * This will be used by the ODE to create system matrices (A = mM + bB + kK). Since it is up to the linear solver
     * to use its matrix format (sparse, full, block, etc.), it must have an API to allow other components to create
     * a compatible matrix that will later be used by this same linear solver.
     *
     * @param rows Number of rows of the matrix.
     * @param cols Number of columns of the matrix.
     * @return A pointer to the newly created matrix.
     *
     * @note Important: The caller of this function (for example, the ODE component) is responsible to free the
     *       memory of this newly created matrix.
     */
    [[nodiscard]]
    virtual BaseMatrix * create_new_matrix(unsigned int rows, unsigned int cols) const = 0;

    /**
     * Creates a new BaseVector of size n.
     *
     * This will be used by the ODE to create system vectors (both the left-hand side solution vector, and the right-
     * hand side force vector). Since it is up to the linear solver to use its vector format (sparse, full, block,
     * etc.), it must have an API to allow other components to create a compatible vector that will later be used by
     * this same linear solver.
     *
     * @param n The number of scalar stored inside this vector.
     * @return A pointer to the newly created vector
     *
     * @note Important: The caller of this function (for example, the ODE component) is responsible to free the
     *       memory of this newly created vector.
     */
    [[nodiscard]]
    virtual BaseVector * create_new_vector(unsigned int n) const = 0;

    /**
     * Solve the linear system A [X] = F
     *
     * @param F The right-hand side (RHS) vector.
     * @param X The left-hand side (LHS) solution vector.
     *
     * @return True when the system has been successfully solved, false otherwise.
     *
     * @note The vectors must be of the virtual type SofaCaribou::Algebra::EigenVector<Vector>
     * @note LinearSolver::factorize must have been called before this method.
     */
    virtual bool solve(const BaseVector * F,
                       BaseVector * X) = 0;

    /**
     * Analyze the pattern of the given matrix.
     *
     * This is usually done just before factorizing the matrix and will often produce a permutation
     * matrix. If the pattern of the matrix doesn't change much between factorizations, it can be
     * beneficial to not call this method too often.
     *
     * @param A The matrix to analyze.
     * @return True if the matrix pattern was successfully analyzed, false otherwise.
     */
    virtual bool analyze_pattern() = 0;

    /**
     * Factorize the given matrix.
     *
     * This should be done before a call to solve.
     *
     * @param A The matrix to factorize.
     * @return True if the matrix was successfully factorized, false otherwise.
     */
    virtual bool factorize() = 0;

    /**
     * Pass the (assembled) system matrix that will have to be solved
     * later on using the analyse_patern(), factorize(), and solve() methods.
     * @param A A pointer to the fully assembled system matrix.
     */
    virtual void set_system_matrix(const BaseMatrix * A) = 0;

    /**
     * Returns true if this solver is an iterative one, false otherwise. In case of iterative solvers,
     * the method squared_residuals() returns the list of squared residual norms of every
     * iterations of the last solve call.
     */
    [[nodiscard]] virtual bool is_iterative() const = 0;

    /**
     * For iterative solvers, returns the list of squared residual norms (||r||^2)
     * of every iterations of the last solve call.
     * For direct solvers, returns an empty set.
     */
    [[nodiscard]] virtual auto squared_residuals() const -> std::vector<FLOATING_POINT_TYPE> = 0;

};

} // namespace SofaCaribou::solver
