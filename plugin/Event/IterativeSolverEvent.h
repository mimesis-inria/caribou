#ifndef CARIBOU_ITERATIVESOLVEREVENT_H
#define CARIBOU_ITERATIVESOLVEREVENT_H

#include <sofa/core/objectmodel/Event.h>

namespace sofa {

namespace caribou {

namespace event {

class IterativeSolverStepBeginEvent : public sofa::core::objectmodel::Event
{
public:

    SOFA_EVENT_H( IterativeSolverStepBeginEvent )

    IterativeSolverStepBeginEvent( );

    const char* getClassName() const override { return "IterativeSolverStepBeginEvent"; }
};

class IterativeSolverStepEndEvent : public sofa::core::objectmodel::Event
{
public:

SOFA_EVENT_H( IterativeSolverStepEndEvent )

    IterativeSolverStepEndEvent( );

    const char* getClassName() const override { return "IterativeSolverStepEndEvent"; }
};

} // namespace event

} // namespace caribou

} // namespace sofa


#endif //CARIBOU_ITERATIVESOLVEREVENT_H
