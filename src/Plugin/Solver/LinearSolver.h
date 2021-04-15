#pragma once

#include <SofaCaribou/config.h>

namespace sofa::defaulttype {
class BaseMatrix;
class BaseVector;
}

namespace SofaCaribou::solver {

/**
 * Base interface for linear solvers. This interface define the generic API that linear solvers in Caribou must
 * provide.
 */
class CARIBOU_API LinearSolver {
public:
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
    virtual sofa::defaulttype::BaseMatrix * create_new_matrix(unsigned int rows, unsigned int cols) const = 0;

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
    virtual sofa::defaulttype::BaseVector * create_new_vector(unsigned int n) const = 0;

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
    virtual bool solve(const sofa::defaulttype::BaseVector * F,
                       sofa::defaulttype::BaseVector * X) const = 0;

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
    virtual bool analyze_pattern(const sofa::defaulttype::BaseMatrix * A) = 0;

    /**
     * Factorize the given matrix.
     *
     * This should be done before a call to solve.
     *
     * @param A The matrix to factorize.
     * @return True if the matrix was successfully factorized, false otherwise.
     */
    virtual bool factorize(const sofa::defaulttype::BaseMatrix * A) = 0;

};

} // namespace SofaCaribou::solver
