#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <Eigen/Core>
#include <Caribou/constants.h>
#include <Caribou/Geometry/Element.h>
#include <Caribou/Geometry/Triangle.h>
#include <Caribou/Geometry/Triangle6.h>
#include <Caribou/Python/Caribou.h>
#include <Caribou/Python/Geometry/Element.h>

namespace py = pybind11;

namespace caribou::geometry::bindings {

template<typename TriangleType>
void declare_triangle(py::module & m, const std::string & name) {
    using BaseTriangleType = typename TriangleType::Base;

    // BaseTriangle
    std::string base_name = "Base" + name;
    declare_element<TriangleType>(m, base_name);
    py::class_<BaseTriangleType, Element<TriangleType>> b (m, base_name.c_str());

    if constexpr(caribou::geometry::traits<TriangleType>::Dimension == 3) {
        b.def("normal", &BaseTriangleType::normal);
    }

    // Triangle
    py::class_<TriangleType, BaseTriangleType> c (m, name.c_str());
    c.def("__str__", [&name](const TriangleType & s) {
        Eigen::IOFormat f(std::numeric_limits<double>::digits10 + 2, 0, ", ", ", ", "[", "]", "[", "]");
        if (caribou::geometry::traits<TriangleType>::Dimension == 1) {
            f.rowPrefix = "";
            f.rowSuffix = "";
        }
        std::stringstream ss;
        ss << name << " : ";
        ss << s.nodes().format(f);
        return ss.str();
    });
}

void create_triangle(pybind11::module & m) {
    declare_triangle<Triangle<_2D>>(m, "Triangle_2D");
    declare_triangle<Triangle6<_2D>>(m, "Triangle6_2D");

    declare_triangle<Triangle<_3D>>(m, "Triangle_3D");
    declare_triangle<Triangle6<_3D>>(m, "Triangle6_3D");

    m.def("Triangle", [](){
        return py::cast(Triangle<_2D>());
    });
    m.def("Triangle6", [](){
        return py::cast(Triangle<_2D>());
    });

    m.def("Triangle", [](const caribou::bindings::Dimension & dim) {
        if (dim == caribou::bindings::Dimension::_2D) {
            return py::cast(Triangle<_2D>());
        } else {
            return py::cast(Triangle<_3D>());
        }
    }, py::arg("dimension"));

    m.def("Triangle6", [](const caribou::bindings::Dimension & dim) {
        if (dim == caribou::bindings::Dimension::_2D) {
            return py::cast(Triangle6<_2D>());
        } else {
            return py::cast(Triangle6<_3D>());
        }
        }, py::arg("dimension"));

    // Creation from 3 nodes in 2D
    m.def("Triangle", [](const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 2> & nodes) {
        return py::cast(Triangle<_2D>(nodes));
    }, py::arg("nodes"));

    m.def("Triangle6", [](const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 2> & nodes) {
        return py::cast(Triangle6<_2D>(Triangle<_2D>(nodes)));
        }, py::arg("nodes"));

    m.def("Triangle", [](const Eigen::Matrix<FLOATING_POINT_TYPE, 2, 1> & n0, const Eigen::Matrix<FLOATING_POINT_TYPE, 2, 1> & n1, const Eigen::Matrix<FLOATING_POINT_TYPE, 2, 1> & n2) {
        return py::cast(Triangle<_2D>(n0, n1, n2));
    }, py::arg("n0"), py::arg("n1"), py::arg("n2"));

    m.def("Triangle6", [](const Eigen::Matrix<FLOATING_POINT_TYPE, 2, 1> & n0, const Eigen::Matrix<FLOATING_POINT_TYPE, 2, 1> & n1, const Eigen::Matrix<FLOATING_POINT_TYPE, 2, 1> & n2) {
        return py::cast(Triangle6<_2D>(n0, n1, n2));
        }, py::arg("n0"), py::arg("n1"), py::arg("n2"));

    // Creation from 3 nodes in 3D
    m.def("Triangle", [](const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 3> & nodes) {
        return py::cast(Triangle<_3D>(nodes));
    }, py::arg("nodes"));

    m.def("Triangle6", [](const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 3> & nodes) {
        return py::cast(Triangle6<_3D>(Triangle<_3D>(nodes)));
        }, py::arg("nodes"));

    m.def("Triangle", [](const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n0, const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n1, const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n2) {
        return py::cast(Triangle<_3D>(n0, n1, n2));
    }, py::arg("n0"), py::arg("n1"), py::arg("n2"));

    m.def("Triangle6", [](const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n0, const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n1, const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n2) {
        return py::cast(Triangle6<_3D>(n0, n1, n2));
        }, py::arg("n0"), py::arg("n1"), py::arg("n2"));

    // Creation from 6 nodes in 2D
    m.def("Triangle6", [](const Eigen::Matrix<FLOATING_POINT_TYPE, 6, 2> & nodes) {
        return py::cast(Triangle6<_2D>(nodes));
    }, py::arg("nodes"));

    m.def("Triangle6", [](const Eigen::Matrix<FLOATING_POINT_TYPE, 2, 1> & n0, const Eigen::Matrix<FLOATING_POINT_TYPE, 2, 1> & n1,
                         const Eigen::Matrix<FLOATING_POINT_TYPE, 2, 1> & n2, const Eigen::Matrix<FLOATING_POINT_TYPE, 2, 1> & n3,
                         const Eigen::Matrix<FLOATING_POINT_TYPE, 2, 1> & n4, const Eigen::Matrix<FLOATING_POINT_TYPE, 2, 1> & n5) {
        return py::cast(Triangle6<_2D>(n0, n1, n2, n3, n4, n5));
    }, py::arg("n0"), py::arg("n1"), py::arg("n2"), py::arg("n3"), py::arg("n4"), py::arg("n5"));

    // Creation from 6 nodes in 3D
    m.def("Triangle6", [](const Eigen::Matrix<FLOATING_POINT_TYPE, 6, 3> & nodes) {
        return py::cast(Triangle6<_3D>(nodes));
    }, py::arg("nodes"));

    m.def("Triangle6", [](const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n0, const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n1,
                         const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n2, const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n3,
                         const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n4, const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n5) {
        return py::cast(Triangle6<_3D>(n0, n1, n2, n3, n4, n5));
    }, py::arg("n0"), py::arg("n1"), py::arg("n2"), py::arg("n3"), py::arg("n4"), py::arg("n5"));

}

} // namespace caribou::geometry::bindings