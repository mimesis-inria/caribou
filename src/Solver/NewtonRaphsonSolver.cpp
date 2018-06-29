#include "NewtonRaphsonSolver.h"

#include "../Event/IterativeSolverEvent.h"

#include <sofa/core/ObjectFactory.h>
#include <sofa/helper/AdvancedTimer.h>
#include <sofa/simulation/MechanicalOperations.h>
#include <sofa/simulation/VectorOperations.h>
#include <sofa/simulation/PropagateEventVisitor.h>

namespace sofa {

namespace caribou {

namespace odesolver {

using core::VecId;
using namespace sofa::defaulttype;
using namespace core::behavior;

NewtonRaphsonSolver::NewtonRaphsonSolver()
    : f_maxit( initData(&f_maxit,(unsigned) 10,"maxIt","Max iterations"))
    , f_corrTolerance( initData(&f_corrTolerance,double(1e-5),"correctionTolerance","tolerance as the norm of correction |du| in each Newton step"))
    , f_resTolerance( initData(&f_resTolerance,double(1e-5),"residualTolerance","tolerance as the norm of residual |f - K(u)| in each Newton step"))
    , f_convergeOnResidual( initData(&f_convergeOnResidual, false,"convergeOnResidual","use the residual as the convergence criterium (stricter)"))

{}

void NewtonRaphsonSolver::solve(const core::ExecParams* params /* PARAMS FIRST */, double /*dt*/, sofa::core::MultiVecCoordId xResult, sofa::core::MultiVecDerivId /*vResult*/) {
    sofa::simulation::common::VectorOperations vop( params, this->getContext() );
    sofa::simulation::common::MechanicalOperations mop( params, this->getContext() );
    sofa::simulation::common::VisitorExecuteFunc executeVisitor(*this->getContext());

    MultiVecCoord x_start(&vop, core::VecCoordId::position() );
    MultiVecCoord x(&vop, xResult /*core::VecCoordId::position()*/ );
    MultiVecDeriv force( &vop, core::VecDerivId::force() );
    dx.realloc( &vop, true );

    // MO vector dx is not allocated by default, it will seg fault if the CG is used (dx is taken by default) with an IdentityMapping
    MultiVecDeriv tempdx(&vop, core::VecDerivId::dx() ); tempdx.realloc( &vop, true, true );

    //////////////////////////////////////////////////////////////////////////////////

    msg_info() << "======= Starting Newton in actual time step " << this->getTime();
    msg_info() << "(maximum of " << f_maxit.getValue() << " iterations)";

    unsigned n_it=0;
    double dx_norm = -1.0, f_norm = -1.0;
    auto context = dynamic_cast<simulation::Node*>(this->getContext());

    sofa::helper::AdvancedTimer::stepBegin("NewtonRaphsonSolver::Solve");

    while (n_it < f_maxit.getValue())
    {

        {
            sofa::helper::AdvancedTimer::stepBegin("step_"+std::to_string(n_it));
            sofa::caribou::event::IterativeSolverStepBeginEvent ev;
            sofa::simulation::PropagateEventVisitor propagator(params, &ev);
            context->execute(propagator);
        }

        // compute addForce, in mapped: addForce + applyJT (vec)
        force.clear();
        mop.computeForce(force);
        mop.projectResponse(force);
        f_norm = sqrt(force.dot(force));

        if (f_convergeOnResidual.getValue() && f_norm <= f_resTolerance.getValue()) {
            msg_info() << "[CONVERGED] The residual's norm |f - K(x0 + dx)| is smaller than the threshold of " << f_resTolerance;
            break;
        }


        // assemble matrix, CG: does nothing
        // LDL non-mapped: addKToMatrix added to system matrix
        // LDL mapped: addKToMatrix not added to the system matrix, needs mapped FF (TODO)
        core::behavior::MultiMatrix<simulation::common::MechanicalOperations> matrix(&mop);
        matrix = MechanicalMatrix::K * -1.0;

        // for CG: calls iteratively addDForce, mapped:  [applyJ, addDForce, applyJt(vec)]+
        // for LDL: solves the system, everything's already assembled
        matrix.solve(dx, force);

        x.eq( x_start, dx, 1 );
        mop.solveConstraint(x, core::ConstraintParams::POS);

        //propagate positions to mapped nodes: taken from AnimateVisitor::processNodeTopDown executed by the animation loop
        // calls apply, applyJ
        sofa::core::MechanicalParams mp;
        sofa::simulation::MechanicalPropagateOnlyPositionAndVelocityVisitor(&mp).execute(this->getContext()); // propagate the changes to mappings below

        dx_norm = sqrt(dx.dot(dx));


        {
            msg_info() << "Iteration #" << n_it << ": |f - K(x0 + dx)| = " << f_norm << " |dx| = " << dx_norm;
            sofa::helper::AdvancedTimer::valSet("residual", f_norm);
            sofa::helper::AdvancedTimer::valSet("correction", dx_norm);
            sofa::caribou::event::IterativeSolverStepEndEvent ev;
            sofa::simulation::PropagateEventVisitor propagator(params, &ev);
            context->execute(propagator);
            sofa::helper::AdvancedTimer::stepEnd("step_"+std::to_string(n_it));
        }

        if (dx_norm <= this->f_corrTolerance.getValue()) {
            msg_info() << "[CONVERGED] The correction's norm |dx| is smaller than the threshold of " << f_corrTolerance;
            break;
        }

        n_it++;

    } // End while (n_it < f_maxit.getValue())

    if (n_it >= f_maxit.getValue()) {
        msg_info() << "[DIVERGED] The number of iterations reached the threshold of " << f_maxit << " iterations";
    }

    sofa::helper::AdvancedTimer::valSet("nb_iterations", n_it+1);
    sofa::helper::AdvancedTimer::valSet("residual", f_norm);
    sofa::helper::AdvancedTimer::valSet("correction", dx_norm);
    sofa::helper::AdvancedTimer::stepEnd("NewtonRaphsonSolver::Solve");
}


SOFA_DECL_CLASS(NewtonRaphsonSolver)

int NewtonRaphsonSolverClass = core::RegisterObject("Newton Raphson Solver")
    .add< NewtonRaphsonSolver >()
;

} // namespace odesolver

} // namespace caribou

} // namespace sofa
