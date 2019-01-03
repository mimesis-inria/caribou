#include <Caribou/Geometry/Python/Node.h>

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>

#include <Caribou/Geometry/Node.h>

namespace py = pybind11;
using namespace caribou::geometry;

namespace caribou {
namespace geometry {
namespace python {

void create_node(py::module & m) {

    py::class_<Node<1>>(m, "Node1D")
            .def(py::init<>())
            .def(py::init<float>())
            .def(py::init<double>())
            .def_property("x",
                          [](const Node<1> & p) { return p.x(); },
                          [](Node<1> & p, const float & x) { p.x() = x; }
            )
            ;

    py::class_<Node<2>>(m, "Node2D")
            .def(py::init<>())
            .def(py::init<float, float>())
            .def(py::init<double, double>())
            .def_property("x",
                          [](const Node<2> &p) { return p.x(); },
                          [](Node<2> &p, const float & x) { p.x() = x; }
            )
            .def_property("y",
                          [](const Node<2> &p) { return p.y(); },
                          [](Node<2> &p, const float & y) { p.y() = y; }
            )
            ;

    py::class_<Node<3>>(m, "Node3D")
            .def(py::init<>())
            .def(py::init<float, float, float>())
            .def(py::init<double, double, double>())
            .def_property("x",
                          [](const Node<3> &p) { return p.x(); },
                          [](Node<3> &p, const float & x) { p.x() = x; }
            )
            .def_property("y",
                          [](const Node<3> &p) { return p.y(); },
                          [](Node<3> &p, const float & y) { p.y() = y; }
            )
            .def_property("z",
                          [](const Node<3> &p) { return p.z(); },
                          [](Node<3> &p, const float & z) { p.z() = z; }
            )
            ;
}

} // namespace python

} // namespace geometry

} // namespace caribou