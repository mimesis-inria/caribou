#include <iostream>
#include "__init__.h"
#include "objetmodel/__init__.h"

namespace caribou {
namespace python {
namespace sofa {
namespace core {


PyObject * create_module()
{
    CP_MODULE_METHODS_BEGIN(core)
    CP_MODULE_METHODS_END

    CP_MODULE_DOC(core, "");

    PyObject * module = CP_CREATE_MODULE(core);
    PyModule_AddObject(module, "objectmodel", caribou::python::sofa::core::objectmodel::create_module());

    return module;
}

} // namespace core
} // namespace sofa
} // namespace python
} // namespace caribou