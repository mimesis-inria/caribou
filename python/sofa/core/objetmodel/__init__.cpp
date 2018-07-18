#include <iostream>
#include "__init__.h"

namespace caribou {
namespace python {
namespace sofa {
namespace core {
namespace objectmodel {


PyObject * create_module()
{
    CP_MODULE_METHODS_BEGIN(objectmodel)
    CP_MODULE_METHODS_END

    CP_MODULE_DOC(objectmodel, "");

    return CP_CREATE_MODULE(objectmodel);
}

} // namespace objectmodel
} // namespace core
} // namespace sofa
} // namespace python
} // namespace caribou