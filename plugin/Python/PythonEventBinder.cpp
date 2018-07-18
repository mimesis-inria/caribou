#include "PythonEventBinder.h"
#include "../Event/IterativeSolverEvent.h"


#include <sofa/simulation/AnimateBeginEvent.h>
#include <sofa/simulation/AnimateEndEvent.h>
#include <sofa/core/ObjectFactory.h>
#include <SofaPython/PythonEnvironment.h>

enum {
    AnimateBeginEvent,
    AnimateEndEvent,
    IterativeSolverStepBeginEvent,
    IterativeSolverStepEndEvent
};

using namespace sofa::caribou::python;

void PythonEventBinder::handleEvent(sofa::core::objectmodel::Event* event)
{
    if (sofa::simulation::AnimateBeginEvent::checkEventType(event)) {
        raise(AnimateBeginEvent);
    } else if (sofa::simulation::AnimateEndEvent::checkEventType(event)) {
        raise(AnimateEndEvent);
    } else if (sofa::caribou::event::IterativeSolverStepBeginEvent::checkEventType(event)) {
        raise(IterativeSolverStepBeginEvent);
    } else if (sofa::caribou::event::IterativeSolverStepEndEvent::checkEventType(event)) {
        raise(IterativeSolverStepEndEvent);
    }
}

void PythonEventBinder::raise(size_t event_index)
{
    if (event_index >= m_callbacks.size())
        return;

    for (const auto & callback : m_callbacks[event_index]) {
        PyObject_CallObject(callback, nullptr);
    }
}

void PythonEventBinder::bind(size_t event_index, PyObject* o)
{
    if (event_index+1 > m_callbacks.size()) {
        m_callbacks.resize(event_index+1);
    }

    m_callbacks[event_index].push_back(o);
}

void PythonEventBinder::bind(std::string event_name, PyObject* o)
{
    if (event_name == "AnimateBeginEvent") {
        bind(AnimateBeginEvent, o);
    } else if (event_name == "AnimateEndEvent") {
        bind(AnimateEndEvent, o);
    } else if (event_name == "IterativeSolverStepBeginEvent") {
        bind(IterativeSolverStepBeginEvent, o);
    } else if (event_name == "IterativeSolverStepEndEvent") {
        bind(IterativeSolverStepEndEvent, o);
    } else {
        msg_error() << "Trying to bind an undefined event named '" << event_name << "'";
    }
}

SOFA_DECL_CLASS(PythonEventBinder)

int PythonEventBinderClass = sofa::core::RegisterObject("Python event binder")
    .add< PythonEventBinder >()
;