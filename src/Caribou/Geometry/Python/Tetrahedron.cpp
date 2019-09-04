#include <Eigen/Core>

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>

#include <Caribou/Geometry/Tetrahedron.h>
#include "Tetrahedron.h"

namespace py = pybind11;

namespace caribou::geometry::python
{

void create_tetrahedrons(pybind11::module & m) {

    py::class_<Tetrahedron<interpolation::Tetrahedron4>> tet4 (m, "Tetrahedron4");
    tet4.def(py::init<>());
    tet4.def(py::init<Eigen::Ref<const Eigen::Matrix<FLOATING_POINT_TYPE, 4, 3, Eigen::RowMajor>> >());
    tet4.def("nodes", &Tetrahedron<interpolation::Tetrahedron4>::nodes);
    tet4.def("T", &Tetrahedron<interpolation::Tetrahedron4>::T);
    tet4.def("jacobian", &Tetrahedron<interpolation::Tetrahedron4>::jacobian);
    tet4.def_property_readonly_static("gauss_nodes", [](py::object /* self */){
        return Eigen::Map<const Eigen::Matrix<FLOATING_POINT_TYPE, Tetrahedron<interpolation::Tetrahedron4>::number_of_gauss_nodes, 3, Eigen::RowMajor>> (&(Tetrahedron<interpolation::Tetrahedron4>::gauss_nodes[0][0]));
    });
    tet4.def_property_readonly_static("gauss_weights", [](py::object /* self */){
        return Eigen::Map<const Eigen::Matrix<FLOATING_POINT_TYPE, Tetrahedron<interpolation::Tetrahedron4>::number_of_gauss_nodes, 1>> (&(Tetrahedron<interpolation::Tetrahedron4>::gauss_weights[0]));
    });

    py::class_<Tetrahedron<interpolation::Tetrahedron10>> tet10 (m, "Tetrahedron10");
    tet10.def(py::init<>());
    tet10.def(py::init<Eigen::Ref<const Eigen::Matrix<FLOATING_POINT_TYPE, 4, 3, Eigen::RowMajor>> >());
    tet10.def(py::init<Eigen::Ref<const Eigen::Matrix<FLOATING_POINT_TYPE, 10, 3, Eigen::RowMajor>> >());
    tet10.def("can_be_converted_to_linear", &Tetrahedron<interpolation::Tetrahedron10>::can_be_converted_to_linear);
    tet10.def("nodes", &Tetrahedron<interpolation::Tetrahedron10>::nodes);
    tet10.def("T", &Tetrahedron<interpolation::Tetrahedron10>::T);
    tet10.def("jacobian", &Tetrahedron<interpolation::Tetrahedron10>::jacobian);
    tet10.def_property_readonly_static("gauss_nodes", [](py::object /* self */){
        return Eigen::Map<const Eigen::Matrix<FLOATING_POINT_TYPE, Tetrahedron<interpolation::Tetrahedron10>::number_of_gauss_nodes, 3, Eigen::RowMajor>> (&(Tetrahedron<interpolation::Tetrahedron10>::gauss_nodes[0][0]));
    });
    tet10.def_property_readonly_static("gauss_weights", [](py::object /* self */){
        return Eigen::Map<const Eigen::Matrix<FLOATING_POINT_TYPE, Tetrahedron<interpolation::Tetrahedron10>::number_of_gauss_nodes, 1>> (&(Tetrahedron<interpolation::Tetrahedron10>::gauss_weights[0]));
    });
}

} // namespace caribou::geometry::python