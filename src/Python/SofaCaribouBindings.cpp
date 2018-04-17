#include "SofaCaribouBindings.h"

#include <SofaPython/PythonFactory.h>
#include "Bindings/PythonEventBinder.h"
#include "PythonEventBinder.h"

using sofa::PythonFactory;

void bindCaribouPythonModule() {
    SP_ADD_CLASS_IN_FACTORY(PythonEventBinder, sofa::caribou::python::PythonEventBinder);
}