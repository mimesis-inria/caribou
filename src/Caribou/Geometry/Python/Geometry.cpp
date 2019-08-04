#include <pybind11/pybind11.h>

#include "Segment.h"

namespace py = pybind11;
using namespace caribou::geometry::python;

PYBIND11_MODULE(CaribouGeometryPython, m) {
    m.doc() = "Geometry module";

    create_segment(m);
}