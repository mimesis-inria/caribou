#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <Eigen/Core>
#include <Caribou/constants.h>
#include <Caribou/Geometry/Segment.h>
#include <Caribou/Python/Caribou.h>
#include <Caribou/Python/Geometry/Element.h>

namespace py = pybind11;

namespace caribou::geometry::bindings {

static auto order_name(const UNSIGNED_INTEGER_TYPE & order) -> std::string {
    if (order == Linear)
        return "Linear";
    else if (order == Quadratic) {
        return "Quadratic";
    }
    return "";
}

template<UNSIGNED_INTEGER_TYPE Dimension, UNSIGNED_INTEGER_TYPE Order>
void declare_segment(py::module & m) {
    std::string name = "Segment" + std::to_string(Dimension) + "D" + order_name(Order);

    using Segment = Segment<Dimension, Order>;
    using BaseSegment = typename Segment::Base;

    // BaseSegment
    std::string base_name = "Base" + name;
    declare_element<Segment>(m, base_name);
    py::class_<BaseSegment, Element<Segment>> (m, base_name.c_str());

    // Segment
    py::class_<Segment, BaseSegment> c (m, name.c_str());
    c.def("__str__", [](const Segment & s) {
        Eigen::IOFormat f(std::numeric_limits<double>::digits10 + 2, 0, ", ", ", ", "[", "]", "[", "]");
        if (Dimension == 1) {
            f.rowPrefix = "";
            f.rowSuffix = "";
        }
        std::stringstream ss;
        ss << "Segment <" << std::to_string(Dimension) << "D, " << order_name(Order) << "> : ";
        ss << s.nodes().format(f);
        return ss.str();
    });
}

void create_segment(pybind11::module & m) {
    declare_segment<_1D, Linear>(m);
    declare_segment<_1D, Quadratic>(m);
    declare_segment<_2D, Linear>(m);
    declare_segment<_2D, Quadratic>(m);
    declare_segment<_3D, Linear>(m);
    declare_segment<_3D, Quadratic>(m);

    m.def("Segment", [](){
        return py::cast(Segment<_1D, Linear>());
    });

    m.def("Segment", [](const caribou::bindings::Dimension & dim) {
        if (dim == caribou::bindings::Dimension::_1D) {
            return py::cast(Segment<_1D, Linear>());
        } else if (dim == caribou::bindings::Dimension::_2D) {
            return py::cast(Segment<_2D, Linear>());
        } else {
            return py::cast(Segment<_3D, Linear>());
        }
    }, py::arg("dimension"));

    m.def("Segment", [](const caribou::bindings::Order & order) {
        if (order == caribou::bindings::Order::Linear) {
            return py::cast(Segment<_1D, Linear>());
        } else {
            return py::cast(Segment<_1D, Quadratic>());
        }
    }, py::arg("order"));

    m.def("Segment", [](const caribou::bindings::Dimension & dim, const caribou::bindings::Order & order) {
        if (dim == caribou::bindings::Dimension::_1D) {
            if (order == caribou::bindings::Order::Linear)
                return py::cast(Segment<_1D, Linear>());
            else
                return py::cast(Segment<_1D, Quadratic>());
        } else if (dim == caribou::bindings::Dimension::_2D) {
            if (order == caribou::bindings::Order::Linear)
                return py::cast(Segment<_2D, Linear>());
            else
                return py::cast(Segment<_2D, Quadratic>());
        } else {
            if (order == caribou::bindings::Order::Linear)
                return py::cast(Segment<_3D, Linear>());
            else
                return py::cast(Segment<_3D, Quadratic>());
        }
    }, py::arg("dimension"), py::arg("order"));

    // Linear creation
    // 1D
    m.def("Segment", [](const Eigen::Matrix<FLOATING_POINT_TYPE, 2, 1> & nodes) {
        return py::cast(Segment<_1D, Linear>(nodes));
    }, py::arg("nodes"));

    m.def("Segment", [](const FLOATING_POINT_TYPE & n0, const FLOATING_POINT_TYPE & n1) {
        Eigen::Matrix<FLOATING_POINT_TYPE, 2, 1> nodes (n0, n1);
        return py::cast(Segment<_1D, Linear>(nodes));
    }, py::arg("n0"), py::arg("n1"));

    // 2D
    m.def("Segment", [](const Eigen::Matrix<FLOATING_POINT_TYPE, 2, 2> & nodes) {
        return py::cast(Segment<_2D, Linear>(nodes));
    }, py::arg("nodes"));

    m.def("Segment", [](const Eigen::Matrix<FLOATING_POINT_TYPE, 2, 1> & n0, const Eigen::Matrix<FLOATING_POINT_TYPE, 2, 1> & n1) {
        return py::cast(Segment<_2D, Linear>(n0, n1));
    }, py::arg("n0"), py::arg("n1"));

    // 3D
    m.def("Segment", [](const Eigen::Matrix<FLOATING_POINT_TYPE, 2, 3> & nodes) {
        return py::cast(Segment<_3D, Linear>(nodes));
    }, py::arg("nodes"));

    m.def("Segment", [](const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n0, const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n1) {
        return py::cast(Segment<_3D, Linear>(n0, n1));
    }, py::arg("n0"), py::arg("n1"));

    // Quadratic creation
    // 1D
    m.def("Segment", [](const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & nodes) {
        return py::cast(Segment<_1D, Quadratic>(nodes));
    }, py::arg("nodes"));

    m.def("Segment", [](const FLOATING_POINT_TYPE & n0, const FLOATING_POINT_TYPE & n1, const FLOATING_POINT_TYPE & n2) {
        Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> nodes (n0, n1, n2);
        return py::cast(Segment<_1D, Quadratic>(nodes));
    }, py::arg("n0"), py::arg("n1"), py::arg("n2"));

    // 2D
    m.def("Segment", [](const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 2> & nodes) {
        return py::cast(Segment<_2D, Quadratic>(nodes));
    }, py::arg("nodes"));

    m.def("Segment", [](const Eigen::Matrix<FLOATING_POINT_TYPE, 2, 1> & n0, const Eigen::Matrix<FLOATING_POINT_TYPE, 2, 1> & n1, const Eigen::Matrix<FLOATING_POINT_TYPE, 2, 1> & n2) {
        return py::cast(Segment<_2D, Quadratic>(n0, n1, n2));
    }, py::arg("n0"), py::arg("n1"), py::arg("n2"));

    // 3D
    m.def("Segment", [](const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 3> & nodes) {
        return py::cast(Segment<_3D, Quadratic>(nodes));
    }, py::arg("nodes"));

    m.def("Segment", [](const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n0, const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n1, const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n2) {
        return py::cast(Segment<_3D, Quadratic>(n0, n1, n2));
    }, py::arg("n0"), py::arg("n1"), py::arg("n2"));


}

}