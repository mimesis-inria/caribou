#include "Base.h"
#include <sofa/core/objectmodel/Base.h>

namespace caribou {
namespace python {
namespace sofa {
namespace core {
namespace objectmodel {

struct PyBase
{
    PyObject_HEAD
    ::sofa::core::objectmodel::Base::SPtr * object;
};

PyTypeObject * create_base_class()
{

}

} // namespace objectmodel
} // namespace core
} // namespace sofa
} // namespace python
} // namespace caribou