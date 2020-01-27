#ifndef SOFACARIBOU_GRAPHCOMPONENTS_SOLVER_CONJUGATEGRADIENT_H
#define SOFACARIBOU_GRAPHCOMPONENTS_SOLVER_CONJUGATEGRADIENT_H

#include <Caribou/config.h>
#include <sofa/core/behavior/LinearSolver.h>
#include <SofaBaseLinearSolver/FullVector.h>
#include <SofaBaseLinearSolver/FullMatrix.h>
#include <SofaBaseLinearSolver/DefaultMultiMatrixAccessor.h>
#include <sofa/helper/OptionsGroup.h>
#include <Eigen/Sparse>
#include <Eigen/Core>

namespace SofaCaribou::GraphComponents::solver {

using sofa::core::objectmodel::Data;
using namespace sofa::component::linearsolver;
using namespace sofa::core::behavior;

/**
 * Implements a Conjugate Gradient method for solving a linear system of the type Ax = b.
 *
 * No actual system is built by this component, meaning that the matrix A, and the vectors
 * x and b are never accumulated into memory. Instead, this solver relies on the vectors stored
 * in the mechanical objects of the current scene context graph. The multiplication of the matrix
 * A and a given vector x is computed by accumulating the result b=Ax from ff.addDForce(x, b) of
 * every forcefields acting on a given mechanical object.
 */
class ConjugateGradientSolver : public LinearSolver {

public:
    SOFA_CLASS(ConjugateGradientSolver, LinearSolver);

    /// Preconditioning methods
    enum class PreconditioningMethod : unsigned int {
        /// No preconditioning, hence the complete matrix won't be built. (default)
                Identity = 0,

        /// Preconditioning based on the incomplete Cholesky factorization.
                IncompleteCholesky = 1
    };

    /**
     * Reset the complete system (A, x and b are cleared). This does absolutely nothing here since the complete
     * system is never built.
     *
     * This method is called by the MultiMatrix::clear() and MultiMatrix::reset() methods.
     */
    void resetSystem() final {}

    /**
     * Set the linear system matrix A = (mM + bB + kK), storing the coefficients m, b and k of
     * the mechanical M,B,K matrices.
     * No actual matrix is built by this method since it is not required by the Conjugate Gradient algorithm.
     * Only the coefficients m, b and k are stored.
     * @param mparams Contains the coefficients m, b and k of the matrices M, B and K
     */
    void setSystemMBKMatrix(const sofa::core::MechanicalParams* mparams) final;

    /**
     * Gives the identifier of the right-hand side vector b. This identifier will be used to find the actual vector
     * in the mechanical objects of the system.
     */
    void setSystemRHVector(sofa::core::MultiVecDerivId b_id) final;

    /**
     * Gives the identifier of the left-hand side vector x. This identifier will be used to find the actual vector
     * in the mechanical objects of the system.
     */
    void setSystemLHVector(sofa::core::MultiVecDerivId x_id) final;

    /**
     * Solves the system by the conjugate gradient method using the coefficients m, b and k; and
     * the vectors x and b.
     */
    void solveSystem() final;

protected:
    /// Constructor
    ConjugateGradientSolver();

    /// INPUTS
    Data<unsigned int> d_maximum_number_of_iterations;
    Data<FLOATING_POINT_TYPE> d_residual_tolerance_threshold;
    Data< sofa::helper::OptionsGroup > d_preconditioning_method;

private:
    ///< The mechanical parameters containing the m, b and k coefficients.
    const sofa::core::MechanicalParams * p_mechanical_params;

    ///< Accessor used to determine the index of each mechanical object matrix and vector in the global system.
    DefaultMultiMatrixAccessor p_accessor;

    ///< The identifier of the b vector
    sofa::core::MultiVecDerivId p_b_id;

    ///< The identifier of the x vector
    sofa::core::MultiVecDerivId p_x_id;

    ///< Global system matrix (only built when a preconditioning method needs it)
    Eigen::SparseMatrix<FLOATING_POINT_TYPE> p_A;

    ///< Global system solution vector (usually filled with an initial guess or the previous solution)
    Eigen::Matrix<FLOATING_POINT_TYPE, Eigen::Dynamic, 1> p_x;

    ///< Global system right-hand side vector
    Eigen::Matrix<FLOATING_POINT_TYPE, Eigen::Dynamic, 1> p_b;
};

} // namespace SofaCaribou::GraphComponents::solver

#endif //SOFACARIBOU_GRAPHCOMPONENTS_SOLVER_CONJUGATEGRADIENT_H
