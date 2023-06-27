#pragma once

#include <SofaCaribou/config.h>
#include <SofaCaribou/Ode/NewtonRaphsonSolver.h>

DISABLE_ALL_WARNINGS_BEGIN
#include <sofa/core/behavior/DefaultMultiMatrixAccessor.h>
#include <sofa/core/objectmodel/Data.h>
DISABLE_ALL_WARNINGS_BEGIN

namespace SofaCaribou::ode {

/**
 * Implementation of a static ODE solver compatible with non-linear materials.
 *
 * We are trying to solve to following
 * \f{eqnarray*}{
 *     \vect{R}(\vect{x}) - \vect{P} = 0
 * \f}
 *
 * Where \f$\vect{R}\f$ is the (possibly non-linear) internal elastic force residual and \f$\vect{P}\f$ is the external
 * force vector (for example, gravitation force or surface traction).
 *
 * Following the <a href="https://en.wikipedia.org/wiki/Newton's_method#Nonlinear_systems_of_equations">Newton-Raphson method</a>,
 * we pose
 *
 * \f{align*}{
 *     \vect{F}(\vect{x}_{n+1}) &= \vect{R}(\vect{x}_{n+1}) - \vect{P}_n \\
 *     \mat{J} = \frac{\partial \vect{F}}{\partial \vect{x}_{n+1}} \bigg\rvert_{\vect{x}_{n+1}^i} &= \mat{K}(\vect{x}_{n+1})
 * \f}
 *
 * where \f$\vect{x}_{n+1}\f$ is the unknown position vector at the \f$n\f$th time step. We then iteratively solve
 *
 * \f{align*}{
 *     \mat{K}(\vect{x}_{n+1}^i) \left [ \Delta \vect{x}_{n+1}^{i+1} \right ] &= - \vect{F}(\vect{x}_{n+1}^i) \\
 *     \vect{x}_{n+1}^{i+1} &= \vect{x}_{n+1}^{i} + \Delta \vect{x}_{n+1}^{i+1}
 * \f}
 *
 */
class StaticODESolver : public NewtonRaphsonSolver {
public:
    SOFA_CLASS(StaticODESolver, NewtonRaphsonSolver);

    template <typename T>
    using Data = sofa::core::objectmodel::Data<T>;

protected:

    /** @see NewtonRaphsonSolver::assemble_rhs_vector */
    void assemble_rhs_vector(const sofa::core::MechanicalParams & mechanical_parameters,
                             const sofa::core::behavior::MultiMatrixAccessor & matrix_accessor,
                             sofa::core::MultiVecDerivId & f_id,
                             SofaCaribou::Algebra::BaseVector * f) final;

    /** @see NewtonRaphsonSolver::assemble_system_matrix */
    void assemble_system_matrix(const sofa::core::MechanicalParams & mechanical_parameters,
                                sofa::core::behavior::DefaultMultiMatrixAccessor & matrix_accessor,
                                SofaCaribou::Algebra::BaseMatrix * A) final;

    /** @see NewtonRaphsonSolver::propagate_solution_increment */
    void propagate_solution_increment(const sofa::core::MechanicalParams & mechanical_parameters,
                                      const sofa::core::behavior::MultiMatrixAccessor & matrix_accessor,
                                      const SofaCaribou::Algebra::BaseVector * dx,
                                      sofa::core::MultiVecCoordId & x_id,
                                      sofa::core::MultiVecDerivId & v_id,
                                      sofa::core::MultiVecDerivId & dx_id) final;
};

} // namespace SofaCaribou::ode
