#include <pybind11/pybind11.h>
namespace py = pybind11;

namespace caribou::geometry::bindings {
void create_segment(pybind11::module & m);
void create_quad(pybind11::module & m);
void create_triangle(pybind11::module & m);
void create_tetrahedron(pybind11::module & m);
void create_hexahedron(pybind11::module & m);
void create_rectangular_hexahedron(pybind11::module & m);
}

PYBIND11_MODULE(Geometry, m) {
    m.doc() = "Geometry module";

    caribou::geometry::bindings::create_quad(m);
    caribou::geometry::bindings::create_segment(m);
    caribou::geometry::bindings::create_triangle(m);
    caribou::geometry::bindings::create_tetrahedron(m);
    caribou::geometry::bindings::create_hexahedron(m);
    caribou::geometry::bindings::create_rectangular_hexahedron(m);
}