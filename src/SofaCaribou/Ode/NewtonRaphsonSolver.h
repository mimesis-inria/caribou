#pragma once

#include <SofaCaribou/config.h>

DISABLE_ALL_WARNINGS_BEGIN
#include <sofa/core/behavior/OdeSolver.h>
#include <sofa/core/behavior/LinearSolver.h>
#include <sofa/core/objectmodel/Data.h>
#include <SofaCaribou/Algebra/BaseMatrixOperations.h>
#include <SofaCaribou/Algebra/BaseVectorOperations.h>
#include <sofa/core/objectmodel/Link.h>
#include <sofa/helper/OptionsGroup.h>
#include <sofa/core/behavior/DefaultMultiMatrixAccessor.h>
DISABLE_ALL_WARNINGS_END

#include <memory>

namespace SofaCaribou::ode {

/**
 * This class implements a generic <a href="https://en.wikipedia.org/wiki/Newton's_method#Nonlinear_systems_of_equations">Newton-Raphson</a> solver for SOFA.
 *
 * Let \f$\vect{F}(\vect{x})\f$ be a non-linear function. This component will try to solve
 *
 * \f{eqnarray*}{
 *     \vect{F}(\vect{x}) = 0
 * \f}
 *
 * using the following iterative method:
 *
 * \f{align*}{
 *     \mat{J} \left [ \Delta \vect{x}_{n+1}^{i+1} \right ] &= - \vect{F}(\vect{x}_{n+1}^i) \\
 *     \vect{x}_{n+1}^{i+1} &= \vect{x}_{n+1}^{i} + \Delta \vect{x}_{n+1}^{i+1}
 * \f}
 *
 * where \f$\vect{x}_{n}^{i}\f$ is the vector \f$\vect{x}\f$ evaluated at the time step \f$n\f$ and Newton iteration \f$i\f$.
 * \f$\mat{J}\f$ is the jacobian of \f$\vect{F}\f$, and is defined as
 *
 * \f{align*}{
 *     \mat{J} = \frac{\partial \vect{F}}{\partial \vect{x}_{n+1}} \bigg\rvert_{\vect{x}_{n+1}^i}
 * \f}
 *
 *
 */
class NewtonRaphsonSolver : public sofa::core::behavior::OdeSolver {
public:
    SOFA_CLASS(NewtonRaphsonSolver, sofa::core::behavior::OdeSolver);

    template <typename T>
    using Data = sofa::core::objectmodel::Data<T>;

    template <typename T>
    using Link = sofa::core::objectmodel::SingleLink<NewtonRaphsonSolver, T, sofa::core::objectmodel::BaseLink::FLAG_STRONGLINK>;

    /**
     * Different strategies to determine when the pattern of the system matrix should be analyzed in order to,
     * for example, compute a permutation matrix before factorizing it.
     */
    enum class PatternAnalysisStrategy : unsigned int {
        NEVER = 0,
        BEGINNING_OF_THE_SIMULATION,
        BEGINNING_OF_THE_TIME_STEP,
        ALWAYS
    };


    NewtonRaphsonSolver();


    void init() override;


    void reset() override;


    void solve (const sofa::core::ExecParams* params, SReal dt, sofa::core::MultiVecCoordId x_id, sofa::core::MultiVecDerivId v_id) override;

    /** List of times (in nanoseconds) that each Newton-Raphson iteration took to compute in the last call to Solve(). */
    auto iteration_times() const -> const std::vector<UNSIGNED_INTEGER_TYPE> & { return p_times; }

    /** The list of squared residual norms (||r||^2) of every newton iterations of the last solve call. */
    auto squared_residuals() const -> const std::vector<FLOATING_POINT_TYPE> & { return p_squared_residuals; }

    /** The initial squared residual (||r0||^2) of the last solve call. */
    auto squared_initial_residual() const -> const FLOATING_POINT_TYPE & { return p_squared_initial_residual; }

    /** Get the current strategy that determine when the pattern of the system matrix should be analyzed. */

    auto pattern_analysis_strategy() const -> PatternAnalysisStrategy;

    /** Set the current strategy that determine when the pattern of the system matrix should be analyzed. */

    void set_pattern_analysis_strategy(const PatternAnalysisStrategy & strategy);

private:

    /**
     * Compute the right-hand side (RHS) of the equation to be solved.
     * It will be called by the Newton-Raphson solver just before starting the Newton iterations.
     * The resulting RHS (usually the internal + external force vector) must be accumulated into the multi-vector
     * identifier returned by "MechanicalParams::f()".
     *
     * This function **must** be overridden by the specialized ODE solver class (for example, the BackwardEuler).
     *
     * @param mechanical_parameters The set of mechanical parameters defined by the Newton-Raphson. It will contains
     * the identifiers of the x, v, f, dx and df multi-vectors, as well as the current dt.
     * @param matrix_accessor The multi-matrix accessor which contains the current mechanical graph.
     * @param f_id The multi-vector identifier of the right-hand side term (forces).
     * @param f A pointer to the right-hand side vector which has already been allocated and initialized to zero.
     *
     * @note This method is responsible to apply any projective constraints (for example, setting entries to zero for
     *       fixed nodes).
     */
    virtual void assemble_rhs_vector(const sofa::core::MechanicalParams & mechanical_parameters,
                                     const sofa::core::behavior::MultiMatrixAccessor & matrix_accessor,
                                     sofa::core::MultiVecDerivId & f_id,
                                     SofaCaribou::Algebra::BaseVector * f) = 0;

    /**
     * Assemble the left-hand side (LHS) system matrix of the equation to be solved. This will be done at the beginning
     * of each Newton step. The system matrix correspond to the jacobian of a function with respect to the unknown variable.
     * Hence, this will usually be constant, unless a non-linear term such as those found in hyper-elastic materials is
     * found.
     *
     * @param mechanical_parameters The set of mechanical parameters defined by the Newton-Raphson. It will contains
     * the identifiers of the x, v, f, dx and df multi-vectors, as well as the current dt.
     * @param matrix_accessor The multi-matrix accessor which contains the current mechanical graph.
     * @param A A pointer to the system matrix which has already been allocated and initialized to zero.
     *
     * @note This method is responsible to apply any projective constraints (for example, setting entries to zero for
     *       fixed nodes).
     *
     * @todo We currently have to pass a non-const reference to the DefaultMultiMatrixAccessor since the const
     *       base class reference would loose the ability to construct multiple matrices (M, B and K) separately.
     *       We should rewrite this once we have our own mechanical graph container instead of this old and useless
     *       matrix accessor.
     */
    virtual void assemble_system_matrix(const sofa::core::MechanicalParams & mechanical_parameters,
                                        sofa::core::behavior::DefaultMultiMatrixAccessor & matrix_accessor,
                                        SofaCaribou::Algebra::BaseMatrix * A) = 0;

    /**
     * Propagate the newly solved increment vector.
     *
     * This method should update the geometry, i.e. x += dx. In addition, it should propagate this increment to mapped
     * mechanical objects through mechanical mappings. In the case of dynamic solvers, this method is also usually
     * responsible to integrate the new time step.
     *
     * @param mechanical_parameters The set of mechanical parameters defined by the Newton-Raphson. It will contains
     * the identifiers of the x, v, f, dx and df multi-vectors, as well as the current dt.
     * @param matrix_accessor The multi-matrix accessor which contains the current mechanical graph.
     * @param dx A pointer to the left-hand side vector which contains the solution of the system A [dx] = F
     * @param x_id The position multi-vector identifier.
     * @param v_id The velocity multi-vector identifier.
     * @param dx_id The increment multi-vector identifier.
     *
     *
     * @note This method is responsible to apply any positional or velocity constraints.
     */
    virtual void propagate_solution_increment(const sofa::core::MechanicalParams & mechanical_parameters,
                                              const sofa::core::behavior::MultiMatrixAccessor & matrix_accessor,
                                              const SofaCaribou::Algebra::BaseVector * dx,
                                              sofa::core::MultiVecCoordId & x_id,
                                              sofa::core::MultiVecDerivId & v_id,
                                              sofa::core::MultiVecDerivId & dx_id) = 0;

protected:
    /** Check that the linked linear solver is not null and that it implements the SofaCaribou::solver::LinearSolver interface */

    bool has_valid_linear_solver () const;

    /// INPUTS
    Data<unsigned> d_newton_iterations;
    Data<double> d_correction_tolerance_threshold;
    Data<double> d_residual_tolerance_threshold;
    Data<double> d_absolute_residual_tolerance_threshold;
    Data<sofa::helper::OptionsGroup> d_pattern_analysis_strategy;

    Link<sofa::core::behavior::LinearSolver> l_linear_solver;

    /// OUTPUTS
    /// Whether or not the last call to solve converged
    Data<bool> d_converged;

    /// Private members

    /// Global system matrix A = mM + bB + kK
    std::unique_ptr<SofaCaribou::Algebra::BaseMatrix> p_A;

    /// Global system LHS vector (the solution)
    std::unique_ptr<SofaCaribou::Algebra::BaseVector> p_DX;

    /// Global system RHS vector (the forces)
    std::unique_ptr<SofaCaribou::Algebra::BaseVector> p_F;

    /// Total displacement since the beginning of the step
    sofa::core::MultiVecDerivId p_U_id;

    /// List of times (in nanoseconds) took to compute each Newton-Raphson iteration
    std::vector<UNSIGNED_INTEGER_TYPE> p_times;

    /// List of squared residual norms (||r||^2) of every newton iterations of the last solve call.
    std::vector<FLOATING_POINT_TYPE> p_squared_residuals;

    /// Initial squared residual (||r0||^2) of the last solve call.
    FLOATING_POINT_TYPE p_squared_initial_residual {};

    /// Either or not the pattern of the system matrix was analyzed at the beginning of the simulation
    bool p_has_already_analyzed_the_pattern = false;
};
}
