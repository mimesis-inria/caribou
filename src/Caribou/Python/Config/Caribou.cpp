#include <pybind11/pybind11.h>
namespace py = pybind11;

#include <Caribou/Python/Caribou.h>

PYBIND11_MODULE(Caribou, m) {
    m.doc() = "Caribou module";

    py::enum_<caribou::bindings::Dimension>(m, "Dimension")
        .value("_1D", caribou::bindings::Dimension::_1D)
        .value("_2D", caribou::bindings::Dimension::_2D)
        .value("_3D", caribou::bindings::Dimension::_3D)
        .export_values();

    py::enum_<caribou::bindings::Order>(m, "Order")
        .value("Linear", caribou::bindings::Order::Linear)
        .value("Quadratic", caribou::bindings::Order::Quadratic)
        .export_values();
}