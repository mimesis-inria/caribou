#include <pybind11/pybind11.h>
namespace py = pybind11;

namespace caribou::geometry::python {
void create_segment(pybind11::module & m);
void create_quad(pybind11::module & m);
void create_triangle(pybind11::module & m);
void create_tetrahedron(pybind11::module & m);
void create_hexahedron(pybind11::module & m);
void create_rectangular_hexahedron(pybind11::module & m);
}

PYBIND11_MODULE(CaribouGeometryPython, m) {
    m.doc() = "Geometry module";

    caribou::geometry::python::create_quad(m);
    caribou::geometry::python::create_segment(m);
    caribou::geometry::python::create_triangle(m);
    caribou::geometry::python::create_tetrahedron(m);
    caribou::geometry::python::create_hexahedron(m);
    caribou::geometry::python::create_rectangular_hexahedron(m);
}