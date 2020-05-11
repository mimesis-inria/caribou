#include <pybind11/pybind11.h>
namespace py = pybind11;

#include <Caribou/Geometry/Python/Segment.h>

PYBIND11_MODULE(CaribouGeometryPython, m) {
    m.doc() = "Geometry module";

    caribou::geometry::python::create_segment(m);
}