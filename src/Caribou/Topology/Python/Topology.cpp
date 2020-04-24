#include <pybind11/pybind11.h>
#include "HashGrid.h"
#include "Grid.h"
#include "IO/IO.h"

namespace py = pybind11;
using namespace caribou::topology::python;
using namespace caribou::topology::io::python;
PYBIND11_MODULE(CaribouTopologyPython, m) {
    m.doc() = "Topology module";
    create_grid(m);
    create_hashgrid(m);
    create_IO(m);
}