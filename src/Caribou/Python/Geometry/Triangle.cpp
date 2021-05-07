#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <Eigen/Core>
#include <Caribou/constants.h>
#include <Caribou/Geometry/Element.h>
#include <Caribou/Geometry/Segment.h>
#include <Caribou/Geometry/Triangle.h>
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
void declare_triangle(py::module & m) {
    std::string name = "Triangle" + std::to_string(Dimension) + "D" + order_name(Order);

    using Triangle = Triangle<Dimension, Order>;
    using BaseTriangle = typename Triangle::Base;

    // BaseTriangle
    std::string base_name = "Base" + name;
    declare_element<Triangle>(m, base_name);
    py::class_<BaseTriangle, Element<Triangle>> b (m, base_name.c_str());

    if constexpr(Dimension == 3) {
        b.def("normal", &BaseTriangle::normal);
    }

    // Triangle
    py::class_<Triangle, BaseTriangle> c (m, name.c_str());
    c.def("__str__", [](const Triangle & s) {
        Eigen::IOFormat f(std::numeric_limits<double>::digits10 + 2, 0, ", ", ", ", "[", "]", "[", "]");
        if (Dimension == 1) {
            f.rowPrefix = "";
            f.rowSuffix = "";
        }
        std::stringstream ss;
        ss << "Triangle <" << std::to_string(Dimension) << "D, " << order_name(Order) << "> : ";
        ss << s.nodes().format(f);
        return ss.str();
    });
}

void create_triangle(pybind11::module & m) {
    declare_triangle<_2D, Linear>(m);
    declare_triangle<_2D, Quadratic>(m);
    declare_triangle<_3D, Linear>(m);
    declare_triangle<_3D, Quadratic>(m);

    m.def("Triangle", [](){
        return py::cast(Triangle<_2D, Linear>());
    });

    m.def("Triangle", [](const caribou::bindings::Dimension & dim) {
        if (dim == caribou::bindings::Dimension::_2D) {
            return py::cast(Triangle<_2D, Linear>());
        } else {
            return py::cast(Triangle<_3D, Linear>());
        }
    }, py::arg("dimension"));

    m.def("Triangle", [](const caribou::bindings::Order & order) {
        if (order == caribou::bindings::Order::Linear) {
            return py::cast(Triangle<_2D, Linear>());
        } else {
            return py::cast(Triangle<_2D, Quadratic>());
        }
    }, py::arg("order"));

    m.def("Triangle", [](const caribou::bindings::Dimension & dim, const caribou::bindings::Order & order) {
        if (dim == caribou::bindings::Dimension::_2D) {
            if (order == caribou::bindings::Order::Linear)
                return py::cast(Triangle<_2D, Linear>());
            else
                return py::cast(Triangle<_2D, Quadratic>());
        } else {
            if (order == caribou::bindings::Order::Linear)
                return py::cast(Triangle<_3D, Linear>());
            else
                return py::cast(Triangle<_3D, Quadratic>());
        }
    }, py::arg("dimension"), py::arg("order"));

    // Linear creation
    // 2D
    m.def("Triangle", [](const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 2> & nodes, const caribou::bindings::Order & order) {
        return (order == caribou::bindings::Order::Linear)
               ? py::cast(Triangle<_2D, Linear>(nodes))
               : py::cast(Triangle<_2D, Quadratic>(Triangle<_2D, Linear>(nodes)));
    }, py::arg("nodes"), py::arg("order") = caribou::bindings::Order::Linear);

    m.def("Triangle", [](const Eigen::Matrix<FLOATING_POINT_TYPE, 2, 1> & n0, const Eigen::Matrix<FLOATING_POINT_TYPE, 2, 1> & n1, const Eigen::Matrix<FLOATING_POINT_TYPE, 2, 1> & n2,
                         const caribou::bindings::Order & order) {
        return (order == caribou::bindings::Order::Linear)
               ? py::cast(Triangle<_2D, Linear>(n0, n1, n2))
               : py::cast(Triangle<_2D, Quadratic>(n0, n1, n2));
    }, py::arg("n0"), py::arg("n1"), py::arg("n2"), py::arg("order") = caribou::bindings::Order::Linear);

    // 3D
    m.def("Triangle", [](const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 3> & nodes, const caribou::bindings::Order & order) {
        return (order == caribou::bindings::Order::Linear)
               ? py::cast(Triangle<_3D, Linear>(nodes))
               : py::cast(Triangle<_3D, Quadratic>(Triangle<_3D, Linear>(nodes)));
    }, py::arg("nodes"), py::arg("order") = caribou::bindings::Order::Linear);

    m.def("Triangle", [](const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n0, const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n1, const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n2,
                         const caribou::bindings::Order & order) {
        return (order == caribou::bindings::Order::Linear)
               ? py::cast(Triangle<_3D, Linear>(n0, n1, n2))
               : py::cast(Triangle<_3D, Quadratic>(n0, n1, n2));
    }, py::arg("n0"), py::arg("n1"), py::arg("n2"), py::arg("order") = caribou::bindings::Order::Linear);

    // Quadratic creation
    // 2D
    m.def("Triangle", [](const Eigen::Matrix<FLOATING_POINT_TYPE, 6, 2> & nodes) {
        return py::cast(Triangle<_2D, Quadratic>(nodes));
    }, py::arg("nodes"));

    m.def("Triangle", [](const Eigen::Matrix<FLOATING_POINT_TYPE, 2, 1> & n0, const Eigen::Matrix<FLOATING_POINT_TYPE, 2, 1> & n1,
                         const Eigen::Matrix<FLOATING_POINT_TYPE, 2, 1> & n2, const Eigen::Matrix<FLOATING_POINT_TYPE, 2, 1> & n3,
                         const Eigen::Matrix<FLOATING_POINT_TYPE, 2, 1> & n4, const Eigen::Matrix<FLOATING_POINT_TYPE, 2, 1> & n5) {
        return py::cast(Triangle<_2D, Quadratic>(n0, n1, n2, n3, n4, n5));
    }, py::arg("n0"), py::arg("n1"), py::arg("n2"), py::arg("n3"), py::arg("n4"), py::arg("n5"));

    // 3D
    m.def("Triangle", [](const Eigen::Matrix<FLOATING_POINT_TYPE, 6, 3> & nodes) {
        return py::cast(Triangle<_3D, Quadratic>(nodes));
    }, py::arg("nodes"));

    m.def("Triangle", [](const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n0, const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n1,
                         const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n2, const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n3,
                         const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n4, const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n5) {
        return py::cast(Triangle<_3D, Quadratic>(n0, n1, n2, n3, n4, n5));
    }, py::arg("n0"), py::arg("n1"), py::arg("n2"), py::arg("n3"), py::arg("n4"), py::arg("n5"));

}

} // namespace caribou::geometry::bindings