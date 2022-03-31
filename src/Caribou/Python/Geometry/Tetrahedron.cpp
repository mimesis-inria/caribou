#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <Eigen/Core>
#include <Caribou/Geometry/Tetrahedron.h>
#include <Caribou/Geometry/Tetrahedron10.h>
#include <Caribou/Python/Geometry/Element.h>

namespace py = pybind11;

namespace caribou::geometry::bindings {
template<typename TetrahedronType>
void declare_tetrahedron(py::module & m, const std::string & name) {
    using BaseTetrahedron = typename TetrahedronType::Base;

    // BaseTetrahedron
    std::string base_name = "Base" + name;
    declare_element<TetrahedronType>(m, base_name);
    py::class_<BaseTetrahedron, Element<TetrahedronType>> (m, base_name.c_str());

    // Tetrahedron
    py::class_<TetrahedronType, BaseTetrahedron> c (m, name.c_str());
    c.def("__str__", [&name](const TetrahedronType & s) {
        Eigen::IOFormat f(std::numeric_limits<double>::digits10 + 2, 0, ", ", ", ", "[", "]", "[", "]");
        std::stringstream ss;
        ss << name << ": ";
        ss << s.nodes().format(f);
        return ss.str();
    });
}

void create_tetrahedron(pybind11::module & m) {
    declare_tetrahedron<Tetrahedron>(m, "TetrahedronT");
    declare_tetrahedron<Tetrahedron10>(m, "Tetrahedron10T");

    m.def("Tetrahedron", [](){
        return py::cast(Tetrahedron());
    });

    m.def("Tetrahedron10", [](){
        return py::cast(Tetrahedron10());
    });

    // Creation from 4 nodes
    m.def("Tetrahedron", [](const Eigen::Matrix<FLOATING_POINT_TYPE, 4, 3> & nodes) {
        return py::cast(Tetrahedron(nodes));
    }, py::arg("nodes"));

    m.def("Tetrahedron10", [](const Eigen::Matrix<FLOATING_POINT_TYPE, 4, 3> & nodes) {
        return py::cast(Tetrahedron10(Tetrahedron(nodes)));
    }, py::arg("nodes"));

    m.def("Tetrahedron", [](const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n0, const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n1,
                            const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n2, const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n3) {
        return py::cast(Tetrahedron(n0, n1, n2, n3));
    }, py::arg("n0"), py::arg("n1"), py::arg("n2"), py::arg("n3"));

    m.def("Tetrahedron10", [](const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n0, const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n1,
            const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n2, const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n3) {
        return py::cast(Tetrahedron10(n0, n1, n2, n3));
        }, py::arg("n0"), py::arg("n1"), py::arg("n2"), py::arg("n3"));

    // Creation from 10 nodes
    m.def("Tetrahedron10", [](const Eigen::Matrix<FLOATING_POINT_TYPE, 10, 3> & nodes) {
        return py::cast(Tetrahedron10(nodes));
    }, py::arg("nodes"));

    m.def("Tetrahedron10", [](const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n0, const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n1,
                            const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n2, const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n3,
                            const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n4, const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n5,
                            const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n6, const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n7,
                            const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n8, const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n9) {
        return py::cast(Tetrahedron10(n0, n1, n2, n3, n4, n5, n6, n7, n8, n9));
    }, py::arg("n0"), py::arg("n1"), py::arg("n2"), py::arg("n3"), py::arg("n4"), py::arg("n5"), py::arg("n6"), py::arg("n7"), py::arg("n8"), py::arg("n9"));

}

} // namespace caribou::geometry::bindings