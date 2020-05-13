#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <Eigen/Core>
#include <Caribou/Constants.h>
#include <Caribou/Geometry/Segment.h>
#include <Caribou/Geometry/Quad.h>
#include <Caribou/Python/Caribou.h>
#include <Caribou/Geometry/Python/Element.h>

namespace py = pybind11;

namespace caribou::geometry::python {
static auto order_name(const UNSIGNED_INTEGER_TYPE & order) -> std::string {
    if (order == Linear)
        return "Linear";
    else if (order == Quadratic) {
        return "Quadratic";
    }
    return "";
}

template<UNSIGNED_INTEGER_TYPE Dimension, UNSIGNED_INTEGER_TYPE Order>
void declare_quad(py::module & m) {
    std::string name = "Quad" + std::to_string(Dimension) + "D" + order_name(Order);

    using Quad = Quad<Dimension, Order>;
    using BaseQuad = typename Quad::Base;

    // BaseQuad
    std::string base_name = "Base" + name;
    declare_element<Quad>(m, base_name);
    py::class_<BaseQuad, Element<Quad>> (m, base_name.c_str());

    // Quad
    py::class_<Quad, BaseQuad> c (m, name.c_str());
    c.def("__str__", [](const Quad & s) {
        Eigen::IOFormat f(std::numeric_limits<double>::digits10 + 2, 0, ", ", ", ", "[", "]", "[", "]");
        if (Dimension == 1) {
            f.rowPrefix = "";
            f.rowSuffix = "";
        }
        std::stringstream ss;
        ss << "Quad <" << std::to_string(Dimension) << "D, " << order_name(Order) << "> : ";
        ss << s.nodes().format(f);
        return ss.str();
    });
}

void create_quad(pybind11::module & m) {
    declare_quad<_2D, Linear>(m);
    declare_quad<_2D, Quadratic>(m);
    declare_quad<_3D, Linear>(m);
    declare_quad<_3D, Quadratic>(m);

    m.def("Quad", [](){
        return py::cast(Quad<_2D, Linear>());
    });

    m.def("Quad", [](const caribou::python::Dimension & dim) {
        if (dim == caribou::python::Dimension::_2D) {
            return py::cast(Quad<_2D, Linear>());
        } else {
            return py::cast(Quad<_3D, Linear>());
        }
    }, py::arg("dimension"));

    m.def("Quad", [](const caribou::python::Order & order) {
        if (order == caribou::python::Order::Linear) {
            return py::cast(Quad<_2D, Linear>());
        } else {
            return py::cast(Quad<_2D, Quadratic>());
        }
    }, py::arg("order"));

    m.def("Quad", [](const caribou::python::Dimension & dim, const caribou::python::Order & order) {
        if (dim == caribou::python::Dimension::_2D) {
            if (order == caribou::python::Order::Linear)
                return py::cast(Quad<_2D, Linear>());
            else
                return py::cast(Quad<_2D, Quadratic>());
        } else {
            if (order == caribou::python::Order::Linear)
                return py::cast(Quad<_3D, Linear>());
            else
                return py::cast(Quad<_3D, Quadratic>());
        }
    }, py::arg("dimension"), py::arg("order"));

    // Linear creation
    // 2D
    m.def("Quad", [](const Eigen::Matrix<FLOATING_POINT_TYPE, 4, 2> & nodes, const caribou::python::Order & order) {
        return (order == caribou::python::Order::Linear)
               ? py::cast(Quad<_2D, Linear>(nodes))
               : py::cast(Quad<_2D, Quadratic>(Quad<_2D, Linear>(nodes)));
    }, py::arg("nodes"), py::arg("order") = caribou::python::Order::Linear);

    m.def("Quad", [](const Eigen::Matrix<FLOATING_POINT_TYPE, 2, 1> & n0, const Eigen::Matrix<FLOATING_POINT_TYPE, 2, 1> & n1,
                     const Eigen::Matrix<FLOATING_POINT_TYPE, 2, 1> & n2, const Eigen::Matrix<FLOATING_POINT_TYPE, 2, 1> & n3,
                     const caribou::python::Order & order) {
        return (order == caribou::python::Order::Linear)
               ? py::cast(Quad<_2D, Linear>(n0, n1, n2, n3))
               : py::cast(Quad<_2D, Quadratic>(n0, n1, n2, n3));
    }, py::arg("n0"), py::arg("n1"), py::arg("n2"), py::arg("n3"), py::arg("order") = caribou::python::Order::Linear);

    // 3D
    m.def("Quad", [](const Eigen::Matrix<FLOATING_POINT_TYPE, 4, 3> & nodes, const caribou::python::Order & order) {
        return (order == caribou::python::Order::Linear)
               ? py::cast(Quad<_3D, Linear>(nodes))
               : py::cast(Quad<_3D, Quadratic>(Quad<_3D, Linear>(nodes)));
    }, py::arg("nodes"), py::arg("order") = caribou::python::Order::Linear);

    m.def("Quad", [](const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n0, const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n1,
                     const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n2, const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n3,
                     const caribou::python::Order & order) {
        return (order == caribou::python::Order::Linear)
               ? py::cast(Quad<_3D, Linear>(n0, n1, n2, n3))
               : py::cast(Quad<_3D, Quadratic>(n0, n1, n2, n3));
    }, py::arg("n0"), py::arg("n1"), py::arg("n2"), py::arg("n3"), py::arg("order") = caribou::python::Order::Linear);

    // Quadratic creation
    // 2D
    m.def("Quad", [](const Eigen::Matrix<FLOATING_POINT_TYPE, 8, 2> & nodes) {
        return py::cast(Quad<_2D, Quadratic>(nodes));
    }, py::arg("nodes"));

    m.def("Quad", [](const Eigen::Matrix<FLOATING_POINT_TYPE, 2, 1> & n0, const Eigen::Matrix<FLOATING_POINT_TYPE, 2, 1> & n1,
                     const Eigen::Matrix<FLOATING_POINT_TYPE, 2, 1> & n2, const Eigen::Matrix<FLOATING_POINT_TYPE, 2, 1> & n3,
                     const Eigen::Matrix<FLOATING_POINT_TYPE, 2, 1> & n4, const Eigen::Matrix<FLOATING_POINT_TYPE, 2, 1> & n5,
                     const Eigen::Matrix<FLOATING_POINT_TYPE, 2, 1> & n6, const Eigen::Matrix<FLOATING_POINT_TYPE, 2, 1> & n7) {
        return py::cast(Quad<_2D, Quadratic>(n0, n1, n2, n3, n4, n5, n6, n7));
    }, py::arg("n0"), py::arg("n1"), py::arg("n2"), py::arg("n3"), py::arg("n4"), py::arg("n5"), py::arg("n6"), py::arg("n7"));

    // 3D
    m.def("Quad", [](const Eigen::Matrix<FLOATING_POINT_TYPE, 8, 3> & nodes) {
        return py::cast(Quad<_3D, Quadratic>(nodes));
    }, py::arg("nodes"));

    m.def("Quad", [](const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n0, const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n1,
                     const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n2, const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n3,
                     const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n4, const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n5,
                     const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n6, const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n7) {
        return py::cast(Quad<_3D, Quadratic>(n0, n1, n2, n3, n4, n5, n6, n7));
    }, py::arg("n0"), py::arg("n1"), py::arg("n2"), py::arg("n3"), py::arg("n4"), py::arg("n5"), py::arg("n6"), py::arg("n7"));

}

} // namespace caribou::geometry::python