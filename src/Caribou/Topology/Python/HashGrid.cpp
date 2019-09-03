#include <pybind11/pybind11.h>

#include <pybind11/eigen.h>
#include <Caribou/Topology/HashGrid.h>
#include <Caribou/Geometry/Tetrahedron.h>
#include "HashGrid.h"

namespace py = pybind11;

namespace caribou::topology::python
{
using namespace geometry;

template <typename Element>
auto grid(py::module & m, const std::string & name) {
    py::class_<HashGrid<Element, py::object>> o(m, name.c_str());
    o.def("add", &HashGrid<Element, py::object>::add);
    o.def("get", &HashGrid<Element, py::object>::get);
    return o;
}

void create_hashgrid(pybind11::module & m) {
    grid<Tetrahedron<interpolation::Tetrahedron4>>(m, "Tetrahedron4HashGrid");
    grid<Tetrahedron<interpolation::Tetrahedron10>>(m, "Tetrahedron10HashGrid");
}

} // namespace caribou::geometry::python