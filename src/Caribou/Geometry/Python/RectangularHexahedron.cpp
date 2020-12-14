#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <Eigen/Core>
#include <Caribou/Constants.h>
#include <Caribou/Geometry/Quad.h>
#include <Caribou/Geometry/RectangularHexahedron.h>
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

template<UNSIGNED_INTEGER_TYPE Order>
void declare_rectangular_hexahedron(py::module & m) {
    std::string name = "RectangularHexahedron" + order_name(Order);

    using RectangularHexahedron = RectangularHexahedron<Order>;
    using BaseRectangularHexahedron = typename RectangularHexahedron::Base;

    // BaseRectangularHexahedron
    std::string base_name = "Base" + name;
    declare_element<RectangularHexahedron>(m, base_name);
    py::class_<BaseRectangularHexahedron, Element<RectangularHexahedron>> (m, base_name.c_str());

    // RectangularHexahedron
    py::class_<RectangularHexahedron, BaseRectangularHexahedron> c (m, name.c_str());
    c.def("__str__", [](const RectangularHexahedron & s) {
        Eigen::IOFormat vec_f(std::numeric_limits<double>::digits10 + 2, 0, ", ", ", ", "", "", "[", "]");
        Eigen::IOFormat mat_f(std::numeric_limits<double>::digits10 + 2, 0, ", ", ", ", "\n\t  [", "]", "[", "\n\t]");
        std::stringstream ss;
        ss << "RectangularHexahedron <" << order_name(Order) << "> : \n";
        ss << "\tCenter point: " << s.center().format(vec_f) << "\n";
        ss << "\tSize: " << s.size().format(vec_f) << "\n";
        ss << "\tRotation : " << s.rotation().format(mat_f) << "";
        return ss.str();
    });
}

void create_rectangular_hexahedron(pybind11::module & m) {
    declare_rectangular_hexahedron<Linear>(m);
    declare_rectangular_hexahedron<Quadratic>(m);

    using WorldCoordinates = RectangularHexahedron<Linear>::WorldCoordinates;
    using Size = RectangularHexahedron<Linear>::Size;
    using Rotation = RectangularHexahedron<Linear>::Rotation;

    m.def("RectangularHexahedron", [](){
        return py::cast(RectangularHexahedron<Linear>());
    });

    m.def("RectangularHexahedron", [](const caribou::python::Order & order) {
        if (order == caribou::python::Order::Linear) {
            return py::cast(RectangularHexahedron<Linear>());
        } else {
            return py::cast(RectangularHexahedron<Quadratic>());
        }
    }, py::arg("order"));

    m.def("RectangularHexahedron", [](const caribou::python::Order & order, const WorldCoordinates & center) {
        if (order == caribou::python::Order::Linear) {
            return py::cast(RectangularHexahedron<Linear>(center));
        } else {
            return py::cast(RectangularHexahedron<Quadratic>(center));
        }
    }, py::arg("order"), py::arg("center"));

    m.def("RectangularHexahedron", [](const caribou::python::Order & order, const WorldCoordinates & center, const Size & H) {
        if (order == caribou::python::Order::Linear) {
            return py::cast(RectangularHexahedron<Linear>(center, H));
        } else {
            return py::cast(RectangularHexahedron<Quadratic>(center, H));
        }
    }, py::arg("order"), py::arg("center"), py::arg("size"));

    m.def("RectangularHexahedron", [](const caribou::python::Order & order, const WorldCoordinates & center, const Size & H, const Rotation & R) {
        if (order == caribou::python::Order::Linear) {
            return py::cast(RectangularHexahedron<Linear>(center, H, R));
        } else {
            return py::cast(RectangularHexahedron<Quadratic>(center, H, R));
        }
    }, py::arg("order"), py::arg("center"), py::arg("size"), py::arg("rotation"));

    // Linear creation
    m.def("RectangularHexahedron", [](const Eigen::Matrix<FLOATING_POINT_TYPE, 8, 3> & nodes, const caribou::python::Order & order) {
        return (order == caribou::python::Order::Linear)
               ? py::cast(RectangularHexahedron<Linear>(nodes))
               : py::cast(RectangularHexahedron<Quadratic>(Hexahedron<Linear>(nodes)));
    }, py::arg("nodes"), py::arg("order") = caribou::python::Order::Linear);

    m.def("RectangularHexahedron", [](const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n0, const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n1,
                                             const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n2, const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n3,
                                             const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n4, const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n5,
                                             const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n6, const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n7,
                                             const caribou::python::Order & order) {
        return (order == caribou::python::Order::Linear)
               ? py::cast(Hexahedron<Linear>(n0, n1, n2, n3, n4, n5, n6, n7))
               : py::cast(Hexahedron<Quadratic>(n0, n1, n2, n3, n4, n5, n6, n7));
    }, py::arg("n0"), py::arg("n1"), py::arg("n2"), py::arg("n3"), py::arg("n4"), py::arg("n5"), py::arg("n6"), py::arg("n7"), py::arg("order") = caribou::python::Order::Linear);

    // Quadratic creation
    m.def("RectangularHexahedron", [](const Eigen::Matrix<FLOATING_POINT_TYPE, 20, 3> & nodes) {
        return py::cast(RectangularHexahedron<Quadratic>(nodes));
    }, py::arg("nodes"));

    m.def("RectangularHexahedron", [](const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n0,  const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n1,
                                             const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n2,  const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n3,
                                             const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n4,  const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n5,
                                             const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n6,  const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n7,
                                             const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n8,  const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n9,
                                             const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n10, const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n11,
                                             const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n12, const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n13,
                                             const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n14, const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n15,
                                             const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n16, const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n17,
                                             const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n18, const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n19) {
        return py::cast(Hexahedron<Quadratic>(n0, n1, n2, n3, n4, n5, n6, n7, n8, n9, n10, n11, n12, n13, n14, n15, n16, n17, n18, n19));
    }, py::arg("n0"), py::arg("n1"), py::arg("n2"), py::arg("n3"), py::arg("n4"), py::arg("n5"), py::arg("n6"), py::arg("n7"), py::arg("n8"), py::arg("n9"),
       py::arg("n10"), py::arg("n11"), py::arg("n12"), py::arg("n13"), py::arg("n14"), py::arg("n15"), py::arg("n16"), py::arg("n17"), py::arg("n18"), py::arg("n19"));

}

} // namespace caribou::geometry::python