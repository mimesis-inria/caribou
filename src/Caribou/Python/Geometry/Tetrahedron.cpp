#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <Eigen/Core>
#include <Caribou/constants.h>
#include <Caribou/Geometry/Triangle.h>
#include <Caribou/Geometry/Tetrahedron.h>
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

template<UNSIGNED_INTEGER_TYPE Order>
void declare_tetrahedron(py::module & m) {
    std::string name = "Tetrahedron" + order_name(Order);

    using Tetrahedron = Tetrahedron<Order>;
    using BaseTetrahedron = typename Tetrahedron::Base;

    // BaseTetrahedron
    std::string base_name = "Base" + name;
    declare_element<Tetrahedron>(m, base_name);
    py::class_<BaseTetrahedron, Element<Tetrahedron>> (m, base_name.c_str());

    // Tetrahedron
    py::class_<Tetrahedron, BaseTetrahedron> c (m, name.c_str());
    c.def("__str__", [](const Tetrahedron & s) {
        Eigen::IOFormat f(std::numeric_limits<double>::digits10 + 2, 0, ", ", ", ", "[", "]", "[", "]");
        std::stringstream ss;
        ss << "Tetrahedron <" << order_name(Order) << "> : ";
        ss << s.nodes().format(f);
        return ss.str();
    });
}

void create_tetrahedron(pybind11::module & m) {
    declare_tetrahedron<Linear>(m);
    declare_tetrahedron<Quadratic>(m);

    m.def("Tetrahedron", [](){
        return py::cast(Tetrahedron<Linear>());
    });

    m.def("Tetrahedron", [](const caribou::bindings::Order & order) {
        if (order == caribou::bindings::Order::Linear) {
            return py::cast(Tetrahedron<Linear>());
        } else {
            return py::cast(Tetrahedron<Quadratic>());
        }
    }, py::arg("order"));

    // Linear creation
    m.def("Tetrahedron", [](const Eigen::Matrix<FLOATING_POINT_TYPE, 4, 3> & nodes, const caribou::bindings::Order & order) {
        return (order == caribou::bindings::Order::Linear)
               ? py::cast(Tetrahedron<Linear>(nodes))
               : py::cast(Tetrahedron<Quadratic>(Tetrahedron<Linear>(nodes)));
    }, py::arg("nodes"), py::arg("order") = caribou::bindings::Order::Linear);

    m.def("Tetrahedron", [](const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n0, const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n1,
                            const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n2, const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n3,
                            const caribou::bindings::Order & order) {
        return (order == caribou::bindings::Order::Linear)
               ? py::cast(Tetrahedron<Linear>(n0, n1, n2, n3))
               : py::cast(Tetrahedron<Quadratic>(n0, n1, n2, n3));
    }, py::arg("n0"), py::arg("n1"), py::arg("n2"), py::arg("n3"), py::arg("order") = caribou::bindings::Order::Linear);

    // Quadratic creation
    m.def("Tetrahedron", [](const Eigen::Matrix<FLOATING_POINT_TYPE, 10, 3> & nodes) {
        return py::cast(Tetrahedron<Quadratic>(nodes));
    }, py::arg("nodes"));

    m.def("Tetrahedron", [](const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n0, const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n1,
                            const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n2, const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n3,
                            const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n4, const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n5,
                            const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n6, const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n7,
                            const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n8, const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n9) {
        return py::cast(Tetrahedron<Quadratic>(n0, n1, n2, n3, n4, n5, n6, n7, n8, n9));
    }, py::arg("n0"), py::arg("n1"), py::arg("n2"), py::arg("n3"), py::arg("n4"), py::arg("n5"), py::arg("n6"), py::arg("n7"), py::arg("n8"), py::arg("n9"));

}

} // namespace caribou::geometry::bindings