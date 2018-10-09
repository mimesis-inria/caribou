#include <pybind11/pybind11.h>

#include <Caribou/Geometry/Python/Point.h>

namespace py = pybind11;
using namespace caribou::geometry::python;

PYBIND11_MODULE(Geometry, m) {
    m.doc() = "Geometry module";

    create_point(m);

}