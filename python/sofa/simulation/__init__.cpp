#include <iostream>
#include "__init__.h"

namespace caribou {
namespace python {
namespace sofa {
namespace simulation {

static PyObject * test(PyObject *self, PyObject * args) {
    std::cout<<"testing 1 2 3"<<std::endl;
    Py_RETURN_NONE;
}

PyObject * create_module()
{
    CP_MODULE_METHODS_BEGIN(simulation)
    CP_MODULE_METHOD(simulation, test)
    CP_MODULE_METHODS_END

    CP_MODULE_DOC(simulation, "The sofa simulation module");

    return CP_CREATE_MODULE(simulation);
}

} // namespace simulation
} // namespace sofa
} // namespace python
} // namespace caribou