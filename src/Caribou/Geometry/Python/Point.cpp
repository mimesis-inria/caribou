#include <Caribou/Geometry/Python/Point.h>

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>

#include <Caribou/Geometry/Point.h>

namespace py = pybind11;
using namespace caribou::geometry;

namespace caribou
{
namespace geometry
{
namespace python
{

void create_point(py::module & m) {

    py::class_<Point1D>(m, "Point1D")
            .def(py::init<>())
            .def(py::init<float>())
            .def(py::init<double>())
            .def_property("x",
                          [](const Point1D &p) { return p.x(); },
                          [](Point1D &p, const float & x) { p.x() = x; }
            )
            ;

    py::class_<Point2D>(m, "Point2D")
            .def(py::init<>())
            .def(py::init<float, float>())
            .def(py::init<double, double>())
            .def_property("x",
                          [](const Point2D &p) { return p.x(); },
                          [](Point2D &p, const float & x) { p.x() = x; }
            )
            .def_property("y",
                          [](const Point2D &p) { return p.y(); },
                          [](Point2D &p, const float & y) { p.y() = y; }
            )
            ;

    py::class_<Point3D>(m, "Point3D")
            .def(py::init<>())
            .def(py::init<float, float, float>())
            .def(py::init<double, double, double>())
            .def_property("x",
                          [](const Point3D &p) { return p.x(); },
                          [](Point3D &p, const float & x) { p.x() = x; }
            )
            .def_property("y",
                          [](const Point3D &p) { return p.y(); },
                          [](Point3D &p, const float & y) { p.y() = y; }
            )
            .def_property("z",
                          [](const Point3D &p) { return p.z(); },
                          [](Point3D &p, const float & z) { p.z() = z; }
            )
            ;
}

} // namespace python

} // namespace geometry

} // namespace caribou