#include <iostream>
#include "__init__.h"

namespace caribou {
namespace python {
namespace sofa {
namespace simulation {
namespace graph {


PyObject * create_module()
{
    CP_MODULE_METHODS_BEGIN(graph)
    CP_MODULE_METHODS_END

    CP_MODULE_DOC(graph, "The sofa simulation graph module");

    return CP_CREATE_MODULE(graph);
}

} // namespace graph
} // namespace simulation
} // namespace sofa
} // namespace python
} // namespace caribou