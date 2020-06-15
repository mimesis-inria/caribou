#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include <Caribou/Topology/Mesh.h>
#include <Caribou/Python/Caribou.h>
#include <Caribou/Topology/Python/Domain.h>

#include <Caribou/Geometry/Hexahedron.h>
#include <Caribou/Geometry/Quad.h>
#include <Caribou/Geometry/Segment.h>
#include <Caribou/Geometry/Triangle.h>
#include <Caribou/Geometry/Tetrahedron.h>

PYBIND11_MAKE_OPAQUE(std::vector<UNSIGNED_INTEGER_TYPE>)
PYBIND11_MAKE_OPAQUE(std::vector<INTEGER_TYPE>)
PYBIND11_MAKE_OPAQUE(std::vector<float>)
PYBIND11_MAKE_OPAQUE(std::vector<double>)

namespace py = pybind11;

namespace caribou::topology::python {

template <UNSIGNED_INTEGER_TYPE Dim, typename MatrixType, typename Element, typename NodeIndex>
void declare_add_domain(py::class_<Mesh<Dim, MatrixType>> & c) {
    c.def("add_domain", [](Mesh<Dim, MatrixType> & mesh, const std::string & domain_name, const Element &, py::EigenDRef<const Eigen::Matrix<NodeIndex, Eigen::Dynamic, geometry::traits<Element>::NumberOfNodesAtCompileTime>> & node_indices) {
        return mesh.template add_domain<Element, NodeIndex>(domain_name, node_indices);
    }, py::arg("domain_name"), py::arg("element_type"), py::arg("node_indices"));
}

template <UNSIGNED_INTEGER_TYPE Dim, typename MatrixType, typename Element>
void declare_add_domain(py::class_<Mesh<Dim, MatrixType>> & c) {
    declare_add_domain<Dim, MatrixType, Element, int>(c);
    declare_add_domain<Dim, MatrixType, Element, long int>(c);
    declare_add_domain<Dim, MatrixType, Element, long int>(c);
    declare_add_domain<Dim, MatrixType, Element, unsigned int>(c);
    declare_add_domain<Dim, MatrixType, Element, unsigned long int>(c);
    declare_add_domain<Dim, MatrixType, Element, unsigned long int>(c);
}

template <UNSIGNED_INTEGER_TYPE Dim, typename MatrixType>
void declare_add_domain(py::class_<Mesh<Dim, MatrixType>> & c) {
    using namespace geometry;
    declare_add_domain<Dim, MatrixType, Segment<Dim, Linear>>(c);
    declare_add_domain<Dim, MatrixType, Segment<Dim, Quadratic>>(c);
    if constexpr (Dim > 1) {
        declare_add_domain<Dim, MatrixType, Triangle<Dim, Linear>>(c);
        declare_add_domain<Dim, MatrixType, Triangle<Dim, Quadratic>>(c);
        declare_add_domain<Dim, MatrixType, Quad<Dim, Linear>>(c);
        declare_add_domain<Dim, MatrixType, Quad<Dim, Quadratic>>(c);
    }
    if constexpr (Dim > 2) {
        declare_add_domain<Dim, MatrixType, Hexahedron<Linear>>(c);
        declare_add_domain<Dim, MatrixType, Hexahedron<Quadratic>>(c);
        declare_add_domain<Dim, MatrixType, Tetrahedron<Linear>>(c);
        declare_add_domain<Dim, MatrixType, Tetrahedron<Quadratic>>(c);
    }
}

template <UNSIGNED_INTEGER_TYPE Dim, typename MatrixType>
void declare_mesh(py::module & m) {
    std::string name = "Mesh"+std::to_string(Dim)+"D"+typeid(MatrixType).name();
    using M = Mesh<Dim, MatrixType>;
    py::class_<M> c(m, name.c_str());

    c.def("dimension", &M::dimension);
    c.def("number_of_domains", &M::number_of_domains);
    c.def("number_of_nodes", &M::number_of_nodes);

    c.def("domains", &M::domains);
    c.def("domain", [](const M & m, const UNSIGNED_INTEGER_TYPE & i) {return m.domain(i);});
    c.def("domain", [](const M & m, const std::string & name) {return m.domain(name);});

    c.def("position", &M::position);
    c.def("positions", [](const M & m, const std::vector<UNSIGNED_INTEGER_TYPE> & indices) {
        return m.positions(indices);
    }, py::arg("indices").noconvert());
    c.def("positions", [](const M & m, const std::vector<INTEGER_TYPE> & indices) {
        return m.positions(indices);
    }, py::arg("indices").noconvert());
    c.def("positions", [](const M & m, const std::vector<int> & indices) {
        return m.positions(indices);
    }, py::arg("indices").noconvert());

    declare_domains(c);
    declare_add_domain(c);

    m.def("Mesh", [](const MatrixType & nodes) {
        return py::cast(M(nodes));
    }, py::arg("nodes"));
}

template <UNSIGNED_INTEGER_TYPE Dim>
void declare_mesh(py::module &m) {
    declare_mesh<Dim, Eigen::Matrix<float, Eigen::Dynamic, Dim, (Dim>1?Eigen::RowMajor:Eigen::ColMajor)>>(m);
    declare_mesh<Dim, Eigen::Matrix<double, Eigen::Dynamic, Dim, (Dim>1?Eigen::RowMajor:Eigen::ColMajor)>>(m);
}

//template <typename Real>
//auto create_mesh_from_numpy_vector(const py::array & a) {
//    Eigen::Index n = a.shape(0);
//    Eigen::Index stride = a.strides(0) / static_cast<ssize_t>(sizeof(Real));
//
//    if (stride == 1) {
//
//    }
////    Eigen::Map<Eigen::Matrix<Real, Eigen::Dynamic, 1>, 0, Eigen::InnerStride
//}
//
//auto create_mesh_from_numpy_vector(const py::array & a)
//{
//    if (a.dtype().is(py::dtype::of<float>())) {
//        return create_mesh_from_numpy_vector<float>(a);
//    } else if (a.dtype().is(py::dtype::of<double>())) {
//        return create_mesh_from_numpy_vector<double>(a);
//    } else {
//        throw std::runtime_error(
//                std::string("Cannot create a mesh from a vector of type ") + std::string(py::str(a.dtype())));
//    }
//}
//
//template<UNSIGNED_INTEGER_TYPE Dim, typename Real>
//auto create_mesh_from_numpy_array(const py::array & a) {
//    Eigen::Index rows = a.shape(0),
//                 cols = a.shape(1);
//
//}
//
//template<UNSIGNED_INTEGER_TYPE Dim>
//auto create_mesh_from_numpy_array(const py::array & a) {
//    if (a.dtype().is(py::dtype::of<float>())) {
//        return create_mesh_from_numpy_array<Dim, float>(a);
//    } else if (a.dtype().is(py::dtype::of<double>())) {
//        return create_mesh_from_numpy_array<Dim, double>(a);
//    } else {
//        throw std::runtime_error(
//                std::string("Cannot create a mesh from an array of type ") + std::string(py::str(a.dtype())));
//    }
//}

void create_mesh(py::module & m) {
    py::bind_vector<std::vector<UNSIGNED_INTEGER_TYPE>>(m, "VectorULong");
    py::bind_vector<std::vector<INTEGER_TYPE>>(m, "VectorLong");
    py::bind_vector<std::vector<int>>(m, "Vectorint");
    py::bind_vector<std::vector<unsigned int>>(m, "Vectoruint");
    py::bind_vector<std::vector<float>>(m, "VectorFloat");
    py::bind_vector<std::vector<double>>(m, "VectorDouble");

    declare_mesh<1>(m);
    declare_mesh<2>(m);
    declare_mesh<3>(m);

//    m.def("Mesh", [](const py::array & a) {
//        const auto dims = a.ndim();
//        if (dims < 1 || dims > 2)
//            throw std::runtime_error("A mesh can be constructed from a vector of 1D coordinates (D=1), or a "
//                                     "NxD matrix for a mesh of N nodes in dimension D (D=1, 2 or 3)");
//        if (dims == 1) {
//            return create_mesh_from_numpy_vector(a);
//        } else {
//            Eigen::Index np_cols = a.shape(1);
//            if (np_cols == 1) {
//                return create_mesh_from_numpy_array<1>(a);
//            } else if (np_cols == 2) {
//                return create_mesh_from_numpy_array<2>(a);
//            } else if (np_cols == 3) {
//                return create_mesh_from_numpy_array<3>(a);
//            } else {
//                throw std::runtime_error("A mesh can be constructed from a vector of 1D coordinates (D=1), or a "
//                                         "NxD matrix for a mesh of N nodes in dimension D (D=2 or 3)");
//            }
//        }
//    });
}

} // namespace caribou::topology::python