#ifndef CARIBOU_PYTHON_SOFA_SIMULATION_GRAPH_INIT_H
#define CARIBOU_PYTHON_SOFA_SIMULATION_GRAPH_INIT_H

#include <Python.h>
#include <python/sofa/macros.h>

CP_DECLARE_MODULE(graph);

namespace caribou {
namespace python {
namespace sofa {
namespace simulation {
namespace graph {

extern PyObject * create_module();

} // namespace graph
} // namespace simulation
} // namespace sofa
} // namespace python
} // namespace caribou

#endif // CARIBOU_PYTHON_SOFA_SIMULATION_GRAPH_INIT_H
