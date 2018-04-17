#include "IterativeSolverEvent.h"

namespace sofa {

namespace caribou {

namespace event {

SOFA_EVENT_CPP( IterativeSolverStepBeginEvent )
SOFA_EVENT_CPP( IterativeSolverStepEndEvent )

IterativeSolverStepBeginEvent::IterativeSolverStepBeginEvent()
{

}

IterativeSolverStepEndEvent::IterativeSolverStepEndEvent()
{

}

} // namespace event

} // namespace caribou

} // namespace sofa