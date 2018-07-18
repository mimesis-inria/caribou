#include "PythonEventBinder.h"
#include "../PythonEventBinder.h"

#include <SofaPython/PythonToSofa.inl>
#include <SofaPython/PythonFactory.h>

using namespace sofa::defaulttype;

static inline sofa::caribou::python::PythonEventBinder * get_binder(PyObject * obj) {
    return sofa::py::unwrap<sofa::caribou::python::PythonEventBinder>(obj);
}

static PyObject * PythonEventBinder_bind (PyObject *self, PyObject * args) {
    sofa::caribou::python::PythonEventBinder * binder = get_binder(self);
    char * name;
    PyObject * callback;

    if (!binder) {
        SP_MESSAGE_ERROR("Failed to bind the PythonEventBinder class.")
        return Py_None;
    }

    if(!PyArg_ParseTuple(args, "sO", &name, &callback)) {
        SP_MESSAGE_ERROR("The bind method needs two arguments: an event 'name' (string) and a python method callback method")
        return Py_None;
    }

    if (!strlen(name)) {
        SP_MESSAGE_ERROR("The bind method needs a valid event name as first argument")
        return Py_None;
    }

    Py_INCREF(callback);
    binder->bind(std::string(name), callback);

    return Py_None;
}

SP_CLASS_METHODS_BEGIN(PythonEventBinder)
SP_CLASS_METHOD(PythonEventBinder, bind)
SP_CLASS_METHODS_END;

SP_CLASS_TYPE_SPTR(PythonEventBinder, sofa::caribou::python::PythonEventBinder, BaseObject)
