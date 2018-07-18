#ifndef CARIBOU_PYTHON_SOFA_CORE_INIT_H
#define CARIBOU_PYTHON_SOFA_CORE_INIT_H

#include <Python.h>
#include <python/sofa/macros.h>

CP_DECLARE_MODULE(core);

namespace caribou {
namespace python {
namespace sofa {
namespace core {

extern PyObject * create_module();

} // namespace core
} // namespace sofa
} // namespace python
} // namespace caribou

#endif // CARIBOU_PYTHON_SOFA_CORE_INIT_H
