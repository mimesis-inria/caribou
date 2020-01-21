#include <pybind11/pybind11.h>
#include "HashGrid.h"
#include "Grid.h"

namespace py = pybind11;
using namespace caribou::topology::python;
PYBIND11_MODULE(CaribouTopologyPython, m) {
    m.doc() = "Topology module";
    create_grid(m);
    create_hashgrid(m);
}