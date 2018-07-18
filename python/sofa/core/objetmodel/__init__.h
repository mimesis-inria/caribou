#ifndef CARIBOU_PYTHON_SOFA_CORE_OBJECTMODEL_INIT_H
#define CARIBOU_PYTHON_SOFA_CORE_OBJECTMODEL_INIT_H

#include <Python.h>
#include <python/sofa/macros.h>

CP_DECLARE_MODULE(objectmodel);

namespace caribou {
namespace python {
namespace sofa {
namespace core {
namespace objectmodel {

extern PyObject * create_module();

} // namespace objectmodel
} // namespace core
} // namespace sofa
} // namespace python
} // namespace caribou

#endif // CARIBOU_PYTHON_SOFA_CORE_OBJECTMODEL_INIT_H
