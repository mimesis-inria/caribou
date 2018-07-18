#include <Python.h>

#include "simulation/__init__.h"
#include "core/__init__.h"

extern "C" void initsofa(void)
{
    if ( !Py_IsInitialized() )
    {
        Py_Initialize();
    }

    PyEval_InitThreads();

    PyObject * libraryModule = Py_InitModule("sofa", nullptr);

    PyModule_AddObject(libraryModule, "simulation", caribou::python::sofa::simulation::create_module());
    PyModule_AddObject(libraryModule, "core", caribou::python::sofa::core::create_module());
}