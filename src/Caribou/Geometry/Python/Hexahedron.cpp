#include <Eigen/Core>

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>

#include <Caribou/Geometry/Hexahedron.h>
#include <Caribou/Geometry/RectangularHexahedron.h>
#include "Hexahedron.h"

namespace py = pybind11;

namespace caribou::geometry::python
{

void create_hexahedrons(pybind11::module & m) {

    using Hexahedron8 = Hexahedron<interpolation::Hexahedron8>;
    py::class_<Hexahedron8> hex8 (m, "Hexahedron8");
    hex8.def(py::init<>());
    hex8.def(py::init<Eigen::Ref<const Eigen::Matrix<FLOATING_POINT_TYPE, 8, 3, Eigen::RowMajor>> >());
    hex8.def("nodes", &Hexahedron8::nodes);
    hex8.def("T", &Hexahedron8::T);
    hex8.def("jacobian", &Hexahedron8::jacobian);
    hex8.def_property_readonly_static("gauss_nodes", [](const py::object & /* self */){
        return Eigen::Map<const Eigen::Matrix<FLOATING_POINT_TYPE, Hexahedron8::number_of_gauss_nodes, 3, Eigen::RowMajor>> (&(Hexahedron8::gauss_nodes[0][0]));
    });
    hex8.def_property_readonly_static("gauss_weights", [](const py::object & /* self */){
        return Eigen::Map<const Eigen::Matrix<FLOATING_POINT_TYPE, Hexahedron8::number_of_gauss_nodes, 1>> (&(Hexahedron8::gauss_weights[0]));
    });
    hex8.def_static("interpolate", [](const typename Hexahedron8::LocalCoordinates & p, Eigen::Ref<const Eigen::Matrix<FLOATING_POINT_TYPE, Hexahedron8::NumberOfNodes, 3, Eigen::RowMajor>> values) {
        return Hexahedron8::interpolate(p, values);
    });

    using RectangularHexahedron8 = RectangularHexahedron<interpolation::Hexahedron8>;
    py::class_<RectangularHexahedron8> rec8 (m, "RectangularHexahedron8");
    rec8.def(py::init<>());
    rec8.def(py::init<Eigen::Ref<const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1>> >(), py::arg("center"));
    rec8.def(py::init<Eigen::Ref<const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1>>, Eigen::Ref<const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1>> >(), py::arg("center"), py::arg("size"));
    rec8.def(py::init<Eigen::Ref<const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1>>, Eigen::Ref<const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1>>, Eigen::Ref<const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 3, Eigen::RowMajor>> >(), py::arg("center"), py::arg("size"), py::arg("rotation"));
    rec8.def(py::init([](Eigen::Ref<const Eigen::Matrix<FLOATING_POINT_TYPE, 8, 3, Eigen::RowMajor>> nodes) {
        Hexahedron8 hex(nodes);

        const auto center = hex.T(Hexahedron8::LocalCoordinates(0,0,0));
        const RectangularHexahedron8 ::Size size (
            (nodes.row(0)-nodes.row(1)).norm(),
            (nodes.row(0)-nodes.row(3)).norm(),
            (nodes.row(0)-nodes.row(4)).norm()
        );

        const auto frame = hex.frame();

        return RectangularHexahedron8(center, size, frame);
    }));
    rec8.def("nodes", &RectangularHexahedron8::nodes);
    rec8.def("T", &RectangularHexahedron8::T);
    rec8.def("Tinv", &RectangularHexahedron8::Tinv);
    rec8.def("contains", &RectangularHexahedron8::contains);
    rec8.def("jacobian", py::overload_cast<>(&RectangularHexahedron8::jacobian, py::const_));
    rec8.def("jacobian", py::overload_cast<const RectangularHexahedron8::LocalCoordinates &>(&RectangularHexahedron8::jacobian, py::const_));
    rec8.def_property_readonly_static("gauss_nodes", [](const py::object & /* self */){
        return Eigen::Map<const Eigen::Matrix<FLOATING_POINT_TYPE, RectangularHexahedron8::number_of_gauss_nodes, 3, Eigen::RowMajor>> (&(RectangularHexahedron8::gauss_nodes[0][0]));
    });
    hex8.def_property_readonly_static("gauss_weights", [](const py::object & /* self */){
        return Eigen::Map<const Eigen::Matrix<FLOATING_POINT_TYPE, RectangularHexahedron8::number_of_gauss_nodes, 1>> (&(RectangularHexahedron8::gauss_weights[0]));
    });
    hex8.def_static("interpolate", [](const typename RectangularHexahedron8::LocalCoordinates & p, Eigen::Ref<const Eigen::Matrix<FLOATING_POINT_TYPE, RectangularHexahedron8::NumberOfNodes, 3, Eigen::RowMajor>> values) {
        return RectangularHexahedron8::interpolate(p, values);
    });
}

} // namespace caribou::geometry::python