#pragma once

#include <Caribou/config.h>
#include <sofa/core/behavior/OdeSolver.h>
#include <sofa/simulation/MechanicalMatrixVisitor.h>
#include <sofa/core/behavior/MultiVec.h>

namespace SofaCaribou::ode {

using sofa::core::objectmodel::Data;

class LegacyStaticODESolver : public sofa::core::behavior::OdeSolver
{
public:
    SOFA_CLASS(LegacyStaticODESolver, sofa::core::behavior::OdeSolver);
    LegacyStaticODESolver();

public:
    void solve (const sofa::core::ExecParams* params /* PARAMS FIRST */, double dt, sofa::core::MultiVecCoordId xResult, sofa::core::MultiVecDerivId vResult) override;

    /**
     * List of times (in nanoseconds) that each Newton-Raphson iteration took to compute in the last call to Solve().
     */
    auto iteration_times() const -> const std::vector<UNSIGNED_INTEGER_TYPE> & {
        return p_times;
    }

    /*!
     * The list of squared residual norms (||r||^2) of every newton iterations of the last solve call.
     */
    auto squared_residuals() const -> const std::vector<FLOATING_POINT_TYPE> & {
        return p_squared_residuals;
    }

    /*!
     * The initial squared residual (||r0||^2) of the last solve call.
     */
    auto squared_initial_residual() const -> const FLOATING_POINT_TYPE & {
        return p_squared_initial_residual;
    }

    /*!
     * The list of squared residual norms (||r||^2) of every iterative linear solver iterations,
     * for each newton iterations of the last solve call.
     */
    auto iterative_linear_solver_squared_residuals() const -> const std::vector<std::vector<FLOATING_POINT_TYPE>> & {
        return p_iterative_linear_solver_squared_residuals;
    }

    /*!
     * List of squared right-hand side norms (||b||^2) of every newton iterations before the call
     * to the solve method of the iterative linear solver.
     */
    auto iterative_linear_solver_squared_rhs_norms() const -> const std::vector<FLOATING_POINT_TYPE> & {
        return p_iterative_linear_solver_squared_rhs_norms;
    }

    /// Given a displacement as computed by the linear system inversion, how much will it affect the velocity
    ///
    /// This method is used to compute the compliance for contact corrections
    /// For Euler methods, it is typically dt.
    double getVelocityIntegrationFactor() const override
    {
        return 1.0; // getContext()->getDt();
    }

    /// Given a displacement as computed by the linear system inversion, how much will it affect the position
    ///
    /// This method is used to compute the compliance for contact corrections
    /// For Euler methods, it is typically dt².
    double getPositionIntegrationFactor() const override
    {
        return getPositionIntegrationFactor(getContext()->getDt());
    }

    virtual double getPositionIntegrationFactor(double dt ) const
    {
        return dt;
    }

    /// Given an input derivative order (0 for position, 1 for velocity, 2 for acceleration),
    /// how much will it affect the output derivative of the given order.
    ///
    /// This method is used to compute the compliance for contact corrections.
    /// For example, a backward-Euler dynamic implicit integrator would use:
    /// Input:      x_t  v_t  a_{t+dt}
    /// x_{t+dt}     1    dt  dt^2
    /// v_{t+dt}     0    1   dt
    ///
    /// If the linear system is expressed on s = a_{t+dt} dt, then the final factors are:
    /// Input:      x_t   v_t    a_t  s
    /// x_{t+dt}     1    dt     0    dt
    /// v_{t+dt}     0    1      0    1
    /// a_{t+dt}     0    0      0    1/dt
    /// The last column is returned by the getSolutionIntegrationFactor method.
    double getIntegrationFactor(int inputDerivative, int outputDerivative) const override
    {
        return getIntegrationFactor(inputDerivative, outputDerivative, getContext()->getDt());
    }

    double getIntegrationFactor(int inputDerivative, int outputDerivative, double dt) const
    {
        double matrix[3][3] =
            {
                { 1, dt, 0},
                { 0, 1, 0},
                { 0, 0, 0}
            };
        if (inputDerivative >= 3 || outputDerivative >= 3)
            return 0;
        else
            return matrix[outputDerivative][inputDerivative];
    }

    /// Given a solution of the linear system,
    /// how much will it affect the output derivative of the given order.
    double getSolutionIntegrationFactor(int outputDerivative) const override
    {
        return getSolutionIntegrationFactor(outputDerivative, getContext()->getDt());
    }

    double getSolutionIntegrationFactor(int outputDerivative, double dt) const
    {
        double vect[3] = { dt, 1, 1/dt};
        if (outputDerivative >= 3)
            return 0;
        else
            return vect[outputDerivative];
    }

protected:

    /// INPUTS
    Data<unsigned> d_newton_iterations;
    Data<double> d_correction_tolerance_threshold;
    Data<double> d_residual_tolerance_threshold;
    Data<bool> d_shoud_diverge_when_residual_is_growing;
    Data<bool> d_warm_start;

    /// OUTPUTS
    ///< Whether or not the last call to solve converged
    Data<bool> d_converged;

private:
    /// Increment at current newton iteration
    sofa::core::behavior::MultiVecDeriv dx;

    /// Total displacement since the beginning of the step
    sofa::core::behavior::MultiVecDeriv U;

    ///< List of times (in nanoseconds) took to compute each Newton-Raphson iteration
    std::vector<UNSIGNED_INTEGER_TYPE> p_times;

    ///< List of squared residual norms (||r||^2) of every newton iterations of the last solve call.
    std::vector<FLOATING_POINT_TYPE> p_squared_residuals;

    ///< Initial squared residual (||r0||^2) of the last solve call.
    FLOATING_POINT_TYPE p_squared_initial_residual;

    ///< List of squared residual norms (||r||^2) of every iterative linear solver iterations,
    ///< for each newton iterations of the last solve call.
    std::vector<std::vector<FLOATING_POINT_TYPE>> p_iterative_linear_solver_squared_residuals;

    ///< List of squared right-hand side norms (||b||^2) of every newton iterations before the call
    ///< to the solve method of the iterative linear solver.
    std::vector<FLOATING_POINT_TYPE> p_iterative_linear_solver_squared_rhs_norms;
};


} // namespace SofaCaribou::ode

