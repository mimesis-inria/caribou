#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <Eigen/Core>
#include <Caribou/Geometry/Quad.h>
#include <Caribou/Geometry/RectangularHexahedron.h>
#include <Caribou/Geometry/RectangularHexahedron20.h>
#include <Caribou/Python/Geometry/Element.h>

namespace py = pybind11;

namespace caribou::geometry::bindings {

template<typename RectangularHexahedronType>
void declare_rectangular_hexahedron(py::module & m, const std::string & name) {
    using BaseRectangularHexahedron = typename RectangularHexahedronType::Base;

    // BaseRectangularHexahedron
    std::string base_name = "Base" + name;
    declare_element<RectangularHexahedronType>(m, base_name);
    py::class_<BaseRectangularHexahedron, Element<RectangularHexahedronType>> (m, base_name.c_str());

    // RectangularHexahedron
    py::class_<RectangularHexahedronType, BaseRectangularHexahedron> c (m, name.c_str());
    c.def("__str__", [&name](const RectangularHexahedronType & s) {
        Eigen::IOFormat vec_f(std::numeric_limits<double>::digits10 + 2, 0, ", ", ", ", "", "", "[", "]");
        Eigen::IOFormat mat_f(std::numeric_limits<double>::digits10 + 2, 0, ", ", ", ", "\n\t  [", "]", "[", "\n\t]");
        std::stringstream ss;
        ss << name << ": \n";
        ss << "\tCenter point: " << s.center().format(vec_f) << "\n";
        ss << "\tSize: " << s.size().format(vec_f) << "\n";
        ss << "\tRotation : " << s.rotation().format(mat_f) << "";
        return ss.str();
    });
}

void create_rectangular_hexahedron(pybind11::module & m) {
    declare_rectangular_hexahedron<RectangularHexahedron>(m, "RectangularHexahedronT");
    declare_rectangular_hexahedron<RectangularHexahedron20>(m, "RectangularHexahedron20T");

    using WorldCoordinates = RectangularHexahedron::WorldCoordinates;
    using Size = RectangularHexahedron::Size;
    using Rotation = RectangularHexahedron::Rotation;

    m.def("RectangularHexahedron", [](){
        return py::cast(RectangularHexahedron());
    });

    m.def("RectangularHexahedron20", [](){
        return py::cast(RectangularHexahedron());
    });

    m.def("RectangularHexahedron", [](const WorldCoordinates & center) {
        return py::cast(RectangularHexahedron(center));
    }, py::arg("center"));

    m.def("RectangularHexahedron20", [](const WorldCoordinates & center) {
        return py::cast(RectangularHexahedron20(center));
    }, py::arg("center"));

    m.def("RectangularHexahedron", [](const WorldCoordinates & center, const Size & H) {
        return py::cast(RectangularHexahedron(center, H));
    }, py::arg("center"), py::arg("size"));

    m.def("RectangularHexahedron20", [](const WorldCoordinates & center, const Size & H) {
        return py::cast(RectangularHexahedron20(center, H));
    }, py::arg("center"), py::arg("size"));

    m.def("RectangularHexahedron", [](const WorldCoordinates & center, const Size & H, const Rotation & R) {
        return py::cast(RectangularHexahedron(center, H, R));
    }, py::arg("center"), py::arg("size"), py::arg("rotation"));

    m.def("RectangularHexahedron20", [](const WorldCoordinates & center, const Size & H, const Rotation & R) {
        return py::cast(RectangularHexahedron20(center, H, R));
    }, py::arg("center"), py::arg("size"), py::arg("rotation"));


    // Creation from 8 nodes
    m.def("RectangularHexahedron", [](const Eigen::Matrix<FLOATING_POINT_TYPE, 8, 3> & nodes) {
        return py::cast(RectangularHexahedron(nodes));
    }, py::arg("nodes"));

    m.def("RectangularHexahedron20", [](const Eigen::Matrix<FLOATING_POINT_TYPE, 8, 3> & nodes) {
        return py::cast(RectangularHexahedron20(Hexahedron(nodes)));
    }, py::arg("nodes"));

    m.def("RectangularHexahedron", [](const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n0, const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n1,
                                             const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n2, const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n3,
                                             const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n4, const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n5,
                                             const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n6, const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n7) {
        return py::cast(RectangularHexahedron(Hexahedron(n0, n1, n2, n3, n4, n5, n6, n7)));
    }, py::arg("n0"), py::arg("n1"), py::arg("n2"), py::arg("n3"), py::arg("n4"), py::arg("n5"), py::arg("n6"), py::arg("n7"));

    m.def("RectangularHexahedron20", [](const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n0, const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n1,
            const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n2, const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n3,
            const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n4, const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n5,
            const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n6, const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n7) {
        return py::cast(RectangularHexahedron20(Hexahedron(n0, n1, n2, n3, n4, n5, n6, n7)));
    }, py::arg("n0"), py::arg("n1"), py::arg("n2"), py::arg("n3"), py::arg("n4"), py::arg("n5"), py::arg("n6"), py::arg("n7"));

    // Creation from 20 nodes
    m.def("RectangularHexahedron20", [](const Eigen::Matrix<FLOATING_POINT_TYPE, 20, 3> & nodes) {
        return py::cast(RectangularHexahedron20(Hexahedron20(nodes)));
    }, py::arg("nodes"));

    m.def("RectangularHexahedron20", [](const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n0,  const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n1,
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