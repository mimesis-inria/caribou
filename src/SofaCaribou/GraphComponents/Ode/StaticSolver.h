#ifndef SOFACARIBOU_GRAPHCOMPONENTS_ODE_STATICSOLVER_H
#define SOFACARIBOU_GRAPHCOMPONENTS_ODE_STATICSOLVER_H

#include <sofa/core/behavior/OdeSolver.h>
#include <sofa/simulation/MechanicalMatrixVisitor.h>
#include <sofa/core/behavior/MultiVec.h>

namespace SofaCaribou {
namespace GraphComponents {
namespace ode {

using sofa::core::objectmodel::Data;

class StaticSolver : public sofa::core::behavior::OdeSolver
{
public:
    SOFA_CLASS(StaticSolver, sofa::core::behavior::OdeSolver);
    StaticSolver();

public:
    void solve (const sofa::core::ExecParams* params /* PARAMS FIRST */, double dt, sofa::core::MultiVecCoordId xResult, sofa::core::MultiVecDerivId vResult) override;

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

    /// the solution vector is stored for warm-start
    sofa::core::behavior::MultiVecDeriv dx;
    Data<unsigned> d_newton_iterations;
    Data<double> d_correction_tolerance_threshold;
    Data<double> d_residual_tolerance_threshold;
    Data<bool> d_should_converge_on_residual_tolerance_threshold;
    Data<bool> d_shoud_diverge_when_residual_is_growing;
};


} // namespace ode
} // namespace GraphComponents
} // namespace SofaCaribou

#endif //SOFACARIBOU_GRAPHCOMPONENTS_ODE_STATICSOLVER_H
