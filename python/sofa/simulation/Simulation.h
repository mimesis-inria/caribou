#ifndef CARIBOU_PYTHON_SOFA_SIMULATION_SIMULATION_H
#define CARIBOU_PYTHON_SOFA_SIMULATION_SIMULATION_H

#include <Python.h>
#include "../macros.h"

//CP_DECLARE_MODULE(simulation);

namespace caribou {
namespace python {
namespace sofa {
namespace simulation {

extern PyTypeObject * create_simulation_class();

} // namespace simulation
} // namespace sofa
} // namespace python
} // namespace caribou


#endif //CARIBOU_PYTHON_SOFA_SIMULATION_SIMULATION_H
