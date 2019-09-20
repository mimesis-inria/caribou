#include <pybind11/pybind11.h>
#include "HashGrid.h"

namespace py = pybind11;
using namespace caribou::topology::python;
PYBIND11_MODULE(CaribouTopologyPython, m) {
    m.doc() = "Topology module";
    create_hashgrid(m);
}