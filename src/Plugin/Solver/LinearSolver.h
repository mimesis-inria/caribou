#pragma once

namespace sofa::defaulttype {
class BaseMatrix;
class BaseVector;
}

namespace SofaCaribou::solver {

/**
 * Base interface for linear solvers. This interface define the generic API that linear solvers in Caribou must
 * provide.
 */
class LinearSolver {
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
     * @param A The assembled system matrix.
     * @param F The right-hand side (RHS) vector.
     * @param X The left-hand side (LHS) solution vector.
     *
     * @return True when the system has been successfully solved, false otherwise.
     *
     * @note The system matrix and vectors are guaranteed to be of the same virtual type as the one returned by
     *       create_new_matrix() and create_new_vector, respectively, since it is these methods that created them
     *       in the first place.
     */
    virtual bool solve(const sofa::defaulttype::BaseMatrix * A,
                       const sofa::defaulttype::BaseVector * F,
                       sofa::defaulttype::BaseVector * X) const = 0;

};

} // namespace SofaCaribou::solver
