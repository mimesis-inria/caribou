#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include <Caribou/Constants.h>
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
    c.def("add_domain", [](Mesh<Dim, MatrixType> & mesh, const std::string & domain_name, const Element &, const Eigen::Matrix<NodeIndex, Eigen::Dynamic, geometry::traits<Element>::NumberOfNodesAtCompileTime> & node_indices) {
        return mesh.template add_domain<Element, NodeIndex>(domain_name, node_indices);
    }, py::arg("domain_name"), py::arg("element_type"), py::arg("node_indices"));

    c.def("add_domain", [](Mesh<Dim, MatrixType> & mesh, const std::string & domain_name, const Element &, const Eigen::Matrix<NodeIndex, Eigen::Dynamic, Eigen::Dynamic> & node_indices) {
        if (geometry::traits<Element>::NumberOfNodesAtCompileTime != caribou::Dynamic and node_indices.cols() != geometry::traits<Element>::NumberOfNodesAtCompileTime) {
            std::ostringstream ss;
            ss << "You tried to create a domain containing elements of " << geometry::traits<Element>::NumberOfNodesAtCompileTime
               << " nodes, but elements having " << node_indices.cols() << " nodes were found instead. The node indices "
               << "matrix should have NxM indices where N is the number of elements, and M is the number of nodes per element.";
            throw py::key_error(ss.str());
        }

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
}

} // namespace caribou::topology::python