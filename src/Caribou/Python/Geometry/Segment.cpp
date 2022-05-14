#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <Eigen/Core>
#include <Caribou/constants.h>
#include <Caribou/Geometry/Segment.h>
#include <Caribou/Geometry/Segment3.h>
#include <Caribou/Python/Caribou.h>
#include <Caribou/Python/Geometry/Element.h>

namespace py = pybind11;

namespace caribou::geometry::bindings {

template<typename SegmentType>
void declare_segment(py::module & m, const std::string & name) {
    using BaseSegmentType = typename SegmentType::Base;

    // BaseSegment
    std::string base_name = "Base" + name;
    declare_element<SegmentType>(m, base_name);
    py::class_<BaseSegmentType, Element<SegmentType>> (m, base_name.c_str());

    // Segment
    py::class_<SegmentType, BaseSegmentType> c (m, name.c_str());
    c.def("__str__", [&name](const SegmentType & s) {
        Eigen::IOFormat f(std::numeric_limits<double>::digits10 + 2, 0, ", ", ", ", "[", "]", "[", "]");
        if constexpr(caribou::geometry::traits<SegmentType>::Dimension == 1) {
            f.rowPrefix = "";
            f.rowSuffix = "";
        }
        std::stringstream ss;
        ss << name << " : ";
        ss << s.nodes().format(f);
        return ss.str();
    });
}

void create_segment(pybind11::module & m) {
    declare_segment<Segment<_1D>>(m, "Segment_1D");
    declare_segment<Segment3<_1D>>(m, "Segment3_1D");

    declare_segment<Segment<_2D>>(m, "Segment_2D");
    declare_segment<Segment3<_2D>>(m, "Segment3_2D");

    declare_segment<Segment<_3D>>(m, "Segment_3D");
    declare_segment<Segment3<_3D>>(m, "Segment3_3D");

    m.def("Segment", [](){
        return py::cast(Segment<_1D>());
    });

    m.def("Segment3", [](){
        return py::cast(Segment3<_1D>());
    });

    m.def("Segment", [](const caribou::bindings::Dimension & dim) {
        if (dim == caribou::bindings::Dimension::_1D) {
            return py::cast(Segment<_1D>());
        } else if (dim == caribou::bindings::Dimension::_2D) {
            return py::cast(Segment<_2D>());
        } else {
            return py::cast(Segment<_3D>());
        }
    }, py::arg("dimension"));

    m.def("Segment3", [](const caribou::bindings::Dimension & dim) {
        if (dim == caribou::bindings::Dimension::_1D) {
            return py::cast(Segment3<_1D>());
        } else if (dim == caribou::bindings::Dimension::_2D) {
            return py::cast(Segment3<_2D>());
        } else {
            return py::cast(Segment3<_3D>());
        }
        }, py::arg("dimension"));

    // Linear creation
    // 1D
    m.def("Segment", [](const Eigen::Matrix<FLOATING_POINT_TYPE, 2, 1> & nodes) {
        return py::cast(Segment<_1D>(nodes));
    }, py::arg("nodes"));

    m.def("Segment", [](const FLOATING_POINT_TYPE & n0, const FLOATING_POINT_TYPE & n1) {
        Eigen::Matrix<FLOATING_POINT_TYPE, 2, 1> nodes (n0, n1);
        return py::cast(Segment<_1D>(nodes));
    }, py::arg("n0"), py::arg("n1"));

    // 2D
    m.def("Segment", [](const Eigen::Matrix<FLOATING_POINT_TYPE, 2, 2> & nodes) {
        return py::cast(Segment<_2D>(nodes));
    }, py::arg("nodes"));

    m.def("Segment", [](const Eigen::Matrix<FLOATING_POINT_TYPE, 2, 1> & n0, const Eigen::Matrix<FLOATING_POINT_TYPE, 2, 1> & n1) {
        return py::cast(Segment<_2D>(n0, n1));
    }, py::arg("n0"), py::arg("n1"));

    // 3D
    m.def("Segment", [](const Eigen::Matrix<FLOATING_POINT_TYPE, 2, 3> & nodes) {
        return py::cast(Segment<_3D>(nodes));
    }, py::arg("nodes"));

    m.def("Segment", [](const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n0, const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n1) {
        return py::cast(Segment<_3D>(n0, n1));
    }, py::arg("n0"), py::arg("n1"));

    // Quadratic creation
    // 1D
    m.def("Segment3", [](const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & nodes) {
        return py::cast(Segment3<_1D>(nodes));
    }, py::arg("nodes"));

    m.def("Segment3", [](const FLOATING_POINT_TYPE & n0, const FLOATING_POINT_TYPE & n1, const FLOATING_POINT_TYPE & n2) {
        Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> nodes (n0, n1, n2);
        return py::cast(Segment3<_1D>(nodes));
    }, py::arg("n0"), py::arg("n1"), py::arg("n2"));

    // 2D
    m.def("Segment3", [](const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 2> & nodes) {
        return py::cast(Segment3<_2D>(nodes));
    }, py::arg("nodes"));

    m.def("Segment3", [](const Eigen::Matrix<FLOATING_POINT_TYPE, 2, 1> & n0, const Eigen::Matrix<FLOATING_POINT_TYPE, 2, 1> & n1, const Eigen::Matrix<FLOATING_POINT_TYPE, 2, 1> & n2) {
        return py::cast(Segment3<_2D>(n0, n1, n2));
    }, py::arg("n0"), py::arg("n1"), py::arg("n2"));

    // 3D
    m.def("Segment3", [](const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 3> & nodes) {
        return py::cast(Segment3<_3D>(nodes));
    }, py::arg("nodes"));

    m.def("Segment3", [](const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n0, const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n1, const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n2) {
        return py::cast(Segment3<_3D>(n0, n1, n2));
    }, py::arg("n0"), py::arg("n1"), py::arg("n2"));


}

}