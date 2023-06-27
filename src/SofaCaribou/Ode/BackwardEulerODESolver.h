#pragma once

#include <SofaCaribou/config.h>
#include <SofaCaribou/Ode/NewtonRaphsonSolver.h>

DISABLE_ALL_WARNINGS_BEGIN
#include <sofa/core/behavior/DefaultMultiMatrixAccessor.h>
#include <sofa/core/objectmodel/Data.h>
DISABLE_ALL_WARNINGS_END

namespace SofaCaribou::ode {

/**
 * Implementation of an implicit backward euler solver compatible with non-linear materials
 *
 * We are trying to solve to following
 * \f{eqnarray*}{
 *     \mat{M} \ddot{\vect{x}} + \mat{C} \dot{\vect{x}} + \vect{R}(\vect{x}) = \vect{P}
 * \f}
 *
 * Where \f$\mat{M}\f$ is the mass matrix, \f$\mat{C}\f$ is the damping matrix, \f$\vect{R}\f$
 * is the (possibly non-linear) internal elastic force residual and \f$\vect{P}\f$ is the external
 * force vector (for example, gravitation force or surface traction).
 *
 * We first transform this second-order differential equation to a first one by introducing two
 * independant variables:
 *
 * \f{eqnarray*}{
 *     \vect{v} &= \dot{\vect{x}} \\
 *     \vect{a} &= \ddot{\vect{x}}
 * \f}
 *
 * Using the <a href="https://en.wikipedia.org/wiki/Backward_Euler_method">backward Euler scheme</a>, we pose the
 * following approximations:
 *
 * \f{align}{
 *     \vect{x}_{n+1} &= \vect{x}_{n} + h \vect{v}_{n+1} \label{eq:backward_euler_x_step} \\
 *     \vect{v}_{n+1} &= \vect{v}_{n} + h \vect{a}_{n+1} \label{eq:backward_euler_v_step}
 * \f}
 *
 * where \f$h\f$ is the delta time between the steps \f$n\f$ and \f$n+1\f$.
 *
 * Substituting \f$(2)\f$ inside \f$(1)\f$ gives
 *
 * \f{eqnarray*}{
 *     \vect{x}_{n+1} &= \vect{x}_{n} + h \left[ \vect{v}_{n} + h \vect{a}_{n+1} \right] \\
 *                    &= \vect{x}_{n} + h \vect{v}_{n} + h^2 \vect{a}_{n+1}
 * \f}
 *
 * And the problem becomes:
 *
 * \f{eqnarray*}{
 *     \mat{M} \vect{a}_{n+1} + \mat{C} \left[ \vect{v}_{n} + h \vect{a}_{n+1} \right] + \vect{R}(\vect{x}_{n} + h \vect{v}_{n} + h^2 \vect{a}_{n+1}) = \vect{P}_n
 * \f}
 * where \f$\vect{a}_{n+1}\f$ is the vector of unknown accelerations.
 *
 * Finally, assuming  \f$\vect{R}\f$ is non-linear in \f$\vect{x}_{n+1}\f$, we iteratively solve for \f$\vect{a}_{n+1}\f$
 * using the
 * <a href="https://en.wikipedia.org/wiki/Newton's_method#Nonlinear_systems_of_equations">Newton-Raphson method</a> and
 * using the previous approximations to back propagate it inside the velocity and position vectors.
 *
 * Let \f$i\f$ be the Newton iteration number for a given time step \f$n\f$, we pose
 *
 * \f{align*}{
 *     \vect{F}(\vect{a}_{n+1}^i) &= \mat{M} \vect{a}_{n+1}^i + \mat{C} \left[ \vect{v}_{n} + h \vect{a}_{n+1}^i \right] + \vect{R}(\vect{x}_{n} + h \vect{v}_{n} + h^2 \vect{a}_{n+1}^i) - \vect{P}_n \\
 *     \mat{J} = \frac{\partial \vect{F}}{\partial \vect{a}_{n+1}} \bigg\rvert_{\vect{a}_{n+1}^i} &= \mat{M} + h \mat{C} + h^2 \mat{K}(\vect{a}_{n+1}^i)
 * \f}
 *
 * where \f$\vect{x}_{n}\f$ and \f$\vect{x}_{n}\f$ are the position and velocity vectors at the beginning of the time
 * step, respectively, and remains constant throughout the Newton iterations.
 *
 * We then solve for \f$\vect{a}_{n+1}^{i+1}\f$ with
 *
 * \f{align*}{
 *     \mat{J} \left [ \Delta \vect{a}_{n+1}^{i+1} \right ] &= - \vect{F}(\vect{a}_{n+1}^i) \\
 *     \vect{a}_{n+1}^{i+1} &= \vect{a}_{n+1}^{i} + \Delta \vect{a}_{n+1}^{i+1}
 * \f}
 *
 * And propagate back the new acceleration using \f$(1)\f$ and \f$(2)\f$.
 *
 * In addition, this component implicitly adds a Rayleigh's damping matrix \f$\mat{C}_r = r_m \mat{M} + r_k \mat{K}(\vect{x}_{n+1})\f$.
 * We therefore have
 *
 * \f{align*}{
 *     \vect{F}(\vect{a}_{n+1}^i) &= \mat{M} \vect{a}_{n+1}^i + \mat{C} \left[ \vect{v}_{n} + h \vect{a}_{n+1}^i \right] + \vect{R}(\vect{x}_{n} + h \vect{v}_{n} + h^2 \vect{a}_{n+1}^i) - \vect{P}_n \\
 *                                &= \mat{M} \vect{a}_{n+1}^i + (\mat{C}_r+\mat{C}) \left[ \vect{v}_{n} + h \vect{a}_{n+1}^i \right] + \vect{R}(\vect{x}_{n} + h \vect{v}_{n} + h^2 \vect{a}_{n+1}^i) - \vect{P}_n \\
 *                                &= \mat{M} \vect{a}_{n+1}^i + (r_m\mat{M}+r_k\mat{K}) \left[ \vect{v}_{n} + h \vect{a}_{n+1}^i \right] + \mat{C} \left[ \vect{v}_{n} + h \vect{a}_{n+1}^i \right] + \vect{R}(\vect{x}_{n} + h \vect{v}_{n} + h^2 \vect{a}_{n+1}^i) - \vect{P}_n \\
 *                                &= \left[ (1 + hr_m)\mat{M} + h\mat{C} + hr_k\mat{K} \right] \vect{a}_{n+1}^i
 *                                 + \left[ r_m\mat{M} + \mat{C} + r_k \mat{K} \right] \vect{v}_n
 *                                 + \left[ \vect{R}(\vect{x}_{n} + h \vect{v}_{n} + h^2 \vect{a}_{n+1}^i) - \vect{P}_n \right] \\
 *     \mat{J} = \frac{\partial \vect{F}}{\partial \vect{a}_{n+1}} \bigg\rvert_{\vect{a}_{n+1}^i} &= (1 + hr_m)\mat{M} + h \mat{C} + h(h+r_k) \mat{K}(\vect{a}_{n+1}^i)
 * \f}
 *
 *
 *
 */
class BackwardEulerODESolver : public NewtonRaphsonSolver {
public:
    SOFA_CLASS(BackwardEulerODESolver, NewtonRaphsonSolver);

    template <typename T>
    using Data = sofa::core::objectmodel::Data<T>;


    BackwardEulerODESolver();


    void solve (const sofa::core::ExecParams* params, SReal dt, sofa::core::MultiVecCoordId x_id, sofa::core::MultiVecDerivId v_id) override;
private:

    /** @see NewtonRaphsonSolver::assemble_rhs_vector */

    void assemble_rhs_vector(const sofa::core::MechanicalParams & mechanical_parameters,
                             const sofa::core::behavior::MultiMatrixAccessor & matrix_accessor,
                             sofa::core::MultiVecDerivId & f_id,
                             SofaCaribou::Algebra::BaseVector * f) final;

    /** @see NewtonRaphsonSolver::assemble_system_matrix */

    void assemble_system_matrix(const sofa::core::MechanicalParams & mechanical_parameters,
                                sofa::core::behavior::DefaultMultiMatrixAccessor & matrix_accessor,
                                SofaCaribou::Algebra::BaseMatrix * A) final;

    /** @see NewtonRaphsonSolver::propagate_position_increment */

    void propagate_solution_increment(const sofa::core::MechanicalParams & mechanical_parameters,
                                      const sofa::core::behavior::MultiMatrixAccessor & matrix_accessor,
                                      const SofaCaribou::Algebra::BaseVector * dx,
                                      sofa::core::MultiVecCoordId & x_id,
                                      sofa::core::MultiVecDerivId & v_id,
                                      sofa::core::MultiVecDerivId & dx_id) final;

    /// INPUTS
    Data<double> d_rayleigh_stiffness;
    Data<double> d_rayleigh_mass;

    /// Private members

    /// Multi-vector identifier of the positions at the beginning of the time step
    sofa::core::MultiVecCoordId p_previous_x_id;

    /// Multi-vector identifier of the velocities at the beginning of the time step
    sofa::core::MultiVecDerivId p_previous_v_id;

    /// Multi-vector identifier of the acceleration at the current newton iteration
    sofa::core::MultiVecDerivId p_a_id;
};

} // namespace SofaCaribou::ode
