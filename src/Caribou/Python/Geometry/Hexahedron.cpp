#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <Eigen/Core>
#include <Caribou/Geometry/Quad.h>
#include <Caribou/Geometry/Hexahedron.h>
#include <Caribou/Geometry/Hexahedron20.h>
#include <Caribou/Python/Geometry/Element.h>

namespace py = pybind11;

namespace caribou::geometry::bindings {

template<typename HexahedronType>
void declare_hexahedron(py::module & m, const std::string & name) {
    using BaseHexahedron = typename HexahedronType::Base;

    // BaseHexahedron
    std::string base_name = "Base" + name;
    declare_element<HexahedronType>(m, base_name);
    py::class_<BaseHexahedron, Element<HexahedronType>> (m, base_name.c_str());

    // Hexahedron
    py::class_<HexahedronType, BaseHexahedron> c (m, name.c_str());
    c.def("__str__", [&name](const Hexahedron & s) {
        Eigen::IOFormat f(std::numeric_limits<double>::digits10 + 2, 0, ", ", ", ", "[", "]", "[", "]");
        std::stringstream ss;
        ss << name << ": ";
        ss << s.nodes().format(f);
        return ss.str();
    });
}

void create_hexahedron(pybind11::module & m) {
    declare_hexahedron<Hexahedron>(m, "HexahedronT");
    declare_hexahedron<Hexahedron20>(m, "Hexahedron20T");

    m.def("Hexahedron", [](){
        return py::cast(Hexahedron());
    });

    m.def("Hexahedron20", [](){
        return py::cast(Hexahedron20());
    });

    // Creation from 8 nodes
    m.def("Hexahedron", [](const Eigen::Matrix<FLOATING_POINT_TYPE, 8, 3> & nodes) {
        return py::cast(Hexahedron(nodes));
    }, py::arg("nodes"));

    m.def("Hexahedron20", [](const Eigen::Matrix<FLOATING_POINT_TYPE, 8, 3> & nodes) {
        return py::cast(Hexahedron20(Hexahedron(nodes)));
    }, py::arg("nodes"));

    m.def("Hexahedron", [](const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n0, const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n1,
                           const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n2, const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n3,
                           const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n4, const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n5,
                           const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n6, const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n7) {
        return py::cast(Hexahedron(n0, n1, n2, n3, n4, n5, n6, n7));
    }, py::arg("n0"), py::arg("n1"), py::arg("n2"), py::arg("n3"), py::arg("n4"), py::arg("n5"), py::arg("n6"), py::arg("n7"));

    m.def("Hexahedron20", [](const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n0, const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n1,
            const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n2, const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n3,
            const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n4, const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n5,
            const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n6, const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n7) {
        return py::cast(Hexahedron20(n0, n1, n2, n3, n4, n5, n6, n7));
    }, py::arg("n0"), py::arg("n1"), py::arg("n2"), py::arg("n3"), py::arg("n4"), py::arg("n5"), py::arg("n6"), py::arg("n7"));

    // Quadratic creation
    m.def("Hexahedron20", [](const Eigen::Matrix<FLOATING_POINT_TYPE, 20, 3> & nodes) {
        return py::cast(Hexahedron20(nodes));
    }, py::arg("nodes"));

    m.def("Hexahedron20", [](const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n0,  const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n1,
                           const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n2,  const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n3,
                           const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n4,  const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n5,
                           const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n6,  const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n7,
                           const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n8,  const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n9,
                           const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n10, const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n11,
                           const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n12, const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n13,
                           const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n14, const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n15,
                           const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n16, const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n17,
                           const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n18, const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n19) {
        return py::cast(Hexahedron20(n0, n1, n2, n3, n4, n5, n6, n7, n8, n9, n10, n11, n12, n13, n14, n15, n16, n17, n18, n19));
    }, py::arg("n0"), py::arg("n1"), py::arg("n2"), py::arg("n3"), py::arg("n4"), py::arg("n5"), py::arg("n6"), py::arg("n7"), py::arg("n8"), py::arg("n9"),
       py::arg("n10"), py::arg("n11"), py::arg("n12"), py::arg("n13"), py::arg("n14"), py::arg("n15"), py::arg("n16"), py::arg("n17"), py::arg("n18"), py::arg("n19"));

}

} // namespace caribou::geometry::bindings