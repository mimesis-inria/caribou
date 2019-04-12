#include <pybind11/pybind11.h>

#include <Caribou/Geometry/Python/Node.h>

namespace py = pybind11;
using namespace caribou::geometry::python;

PYBIND11_MODULE(CaribouGeometryPython, m) {
    m.doc() = "Geometry module";

    create_node(m);
}