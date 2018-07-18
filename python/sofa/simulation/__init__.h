#ifndef CARIBOU_PYTHON_SOFA_SIMULATION_INIT_H
#define CARIBOU_PYTHON_SOFA_SIMULATION_INIT_H

#include <Python.h>
#include "../macros.h"

CP_DECLARE_MODULE(simulation);

namespace caribou {
namespace python {
namespace sofa {
namespace simulation {

extern PyObject * create_module();

} // namespace simulation
} // namespace sofa
} // namespace python
} // namespace caribou

#endif //CARIBOU_PYTHON_SOFA_SIMULATION_INIT_H
