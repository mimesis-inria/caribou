#include <pybind11/pybind11.h>

#include "Segment.h"
#include "Tetrahedron.h"

namespace py = pybind11;
using namespace caribou::geometry::python;

PYBIND11_MODULE(CaribouGeometryPython, m) {
    m.doc() = "Geometry module";

    create_segments(m);
    create_tetrahedrons(m);
}