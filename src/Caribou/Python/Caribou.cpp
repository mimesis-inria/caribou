#include <pybind11/pybind11.h>
namespace py = pybind11;

#include <Caribou/Python/Caribou.h>

PYBIND11_MODULE(CaribouPython, m) {
    m.doc() = "Caribou module";

    py::enum_<caribou::python::Dimension>(m, "Dimension")
        .value("_1D", caribou::python::Dimension::_1D)
        .value("_2D", caribou::python::Dimension::_2D)
        .value("_3D", caribou::python::Dimension::_3D)
        .export_values();

    py::enum_<caribou::python::Order>(m, "Order")
        .value("Linear", caribou::python::Order::Linear)
        .value("Quadratic", caribou::python::Order::Quadratic)
        .export_values();
}