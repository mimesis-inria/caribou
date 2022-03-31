#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <Eigen/Core>
#include <Caribou/constants.h>
#include <Caribou/Geometry/Segment.h>
#include <Caribou/Geometry/Quad.h>
#include <Caribou/Geometry/Quad8.h>
#include <Caribou/Python/Caribou.h>
#include <Caribou/Python/Geometry/Element.h>

namespace py = pybind11;

namespace caribou::geometry::bindings {
template<typename QuadType>
void declare_quad(py::module & m, const std::string & name) {
    using BaseQuad = typename QuadType::Base;

    // BaseQuad
    std::string base_name = "Base" + name;
    declare_element<QuadType>(m, base_name);
    py::class_<BaseQuad, Element<QuadType>> (m, base_name.c_str());

    // Quad
    py::class_<QuadType, BaseQuad> c (m, name.c_str());
    c.def("__str__", [&name](const QuadType & s) {
        Eigen::IOFormat f(std::numeric_limits<double>::digits10 + 2, 0, ", ", ", ", "[", "]", "[", "]");
        if (caribou::geometry::traits<QuadType>::Dimension == 1) {
            f.rowPrefix = "";
            f.rowSuffix = "";
        }
        std::stringstream ss;
        ss << name << ": ";
        ss << s.nodes().format(f);
        return ss.str();
    });
}

void create_quad(pybind11::module & m) {
    declare_quad<Quad<_2D>>(m, "Quad_2D");
    declare_quad<Quad<_3D>>(m, "Quad_3D");
    declare_quad<Quad8<_2D>>(m, "Quad8_2D");
    declare_quad<Quad8<_3D>>(m, "Quad8_3D");

    m.def("Quad", [](){
        return py::cast(Quad<_2D>());
    });

    m.def("Quad", [](const caribou::bindings::Dimension & dim) {
        if (dim == caribou::bindings::Dimension::_2D) {
            return py::cast(Quad<_2D>());
        } else {
            return py::cast(Quad<_3D>());
        }
    }, py::arg("dimension"));

    m.def("Quad8", [](const caribou::bindings::Dimension & dim) {
        if (dim == caribou::bindings::Dimension::_2D) {
            return py::cast(Quad8<_2D>());
        } else {
            return py::cast(Quad8<_3D>());
        }
        }, py::arg("dimension"));

    // Creation from 4 nodes in 2D
    m.def("Quad", [](const Eigen::Matrix<FLOATING_POINT_TYPE, 4, 2> & nodes) {
        return py::cast(Quad<_2D>(nodes));
    }, py::arg("nodes"));

    m.def("Quad8", [](const Eigen::Matrix<FLOATING_POINT_TYPE, 4, 2> & nodes) {
        return py::cast(Quad8<_2D>(Quad<_2D>(nodes)));
        }, py::arg("nodes"));

    m.def("Quad", [](const Eigen::Matrix<FLOATING_POINT_TYPE, 2, 1> & n0, const Eigen::Matrix<FLOATING_POINT_TYPE, 2, 1> & n1,
                     const Eigen::Matrix<FLOATING_POINT_TYPE, 2, 1> & n2, const Eigen::Matrix<FLOATING_POINT_TYPE, 2, 1> & n3) {
        return py::cast(Quad<_2D>(n0, n1, n2, n3));
    }, py::arg("n0"), py::arg("n1"), py::arg("n2"), py::arg("n3"));

    m.def("Quad", [](const Eigen::Matrix<FLOATING_POINT_TYPE, 2, 1> & n0, const Eigen::Matrix<FLOATING_POINT_TYPE, 2, 1> & n1,
            const Eigen::Matrix<FLOATING_POINT_TYPE, 2, 1> & n2, const Eigen::Matrix<FLOATING_POINT_TYPE, 2, 1> & n3) {
        return py::cast(Quad8<_2D>(n0, n1, n2, n3));
        }, py::arg("n0"), py::arg("n1"), py::arg("n2"), py::arg("n3"));

    // Creation from 4 nodes in 3D
    m.def("Quad", [](const Eigen::Matrix<FLOATING_POINT_TYPE, 4, 3> & nodes) {
        return py::cast(Quad<_3D>(nodes));
        }, py::arg("nodes"));

    m.def("Quad8", [](const Eigen::Matrix<FLOATING_POINT_TYPE, 4, 3> & nodes) {
        return py::cast(Quad8<_3D>(Quad<_3D>(nodes)));
        }, py::arg("nodes"));

    m.def("Quad", [](const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n0, const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n1,
            const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n2, const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n3) {
        return py::cast(Quad<_3D>(n0, n1, n2, n3));
        }, py::arg("n0"), py::arg("n1"), py::arg("n2"), py::arg("n3"));

    m.def("Quad8", [](const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n0, const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n1,
            const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n2, const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n3) {
        return py::cast(Quad8<_3D>(n0, n1, n2, n3));
        }, py::arg("n0"), py::arg("n1"), py::arg("n2"), py::arg("n3"));

    // Creation from 8 nodes in 2D
    m.def("Quad8", [](const Eigen::Matrix<FLOATING_POINT_TYPE, 8, 2> & nodes) {
        return py::cast(Quad8<_2D>(nodes));
    }, py::arg("nodes"));

    m.def("Quad8", [](const Eigen::Matrix<FLOATING_POINT_TYPE, 2, 1> & n0, const Eigen::Matrix<FLOATING_POINT_TYPE, 2, 1> & n1,
                     const Eigen::Matrix<FLOATING_POINT_TYPE, 2, 1> & n2, const Eigen::Matrix<FLOATING_POINT_TYPE, 2, 1> & n3,
                     const Eigen::Matrix<FLOATING_POINT_TYPE, 2, 1> & n4, const Eigen::Matrix<FLOATING_POINT_TYPE, 2, 1> & n5,
                     const Eigen::Matrix<FLOATING_POINT_TYPE, 2, 1> & n6, const Eigen::Matrix<FLOATING_POINT_TYPE, 2, 1> & n7) {
        return py::cast(Quad8<_2D>(n0, n1, n2, n3, n4, n5, n6, n7));
    }, py::arg("n0"), py::arg("n1"), py::arg("n2"), py::arg("n3"), py::arg("n4"), py::arg("n5"), py::arg("n6"), py::arg("n7"));

    // 3D
    m.def("Quad8", [](const Eigen::Matrix<FLOATING_POINT_TYPE, 8, 3> & nodes) {
        return py::cast(Quad8<_3D>(nodes));
    }, py::arg("nodes"));

    m.def("Quad8", [](const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n0, const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n1,
                     const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n2, const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n3,
                     const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n4, const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n5,
                     const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n6, const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n7) {
        return py::cast(Quad8<_3D>(n0, n1, n2, n3, n4, n5, n6, n7));
    }, py::arg("n0"), py::arg("n1"), py::arg("n2"), py::arg("n3"), py::arg("n4"), py::arg("n5"), py::arg("n6"), py::arg("n7"));

}

} // namespace caribou::geometry::bindings