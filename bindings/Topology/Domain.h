#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

#include <Caribou/Geometry/Hexahedron.h>
#include <Caribou/Geometry/Quad.h>
#include <Caribou/Geometry/Segment.h>
#include <Caribou/Geometry/Triangle.h>
#include <Caribou/Geometry/Tetrahedron.h>

#include <Caribou/Topology/Mesh.h>
#include <Caribou/Topology/Domain.h>

#include <Caribou/Bindings/Caribou.h>
#include <Caribou/Bindings/Topology/BarycentricContainer.h>

namespace py = pybind11;

namespace caribou::topology::bindings {
template <typename Mesh, typename Element, typename NodeIndex>
void declare_domain(py::class_<Mesh> & m, const std::string & name) {
    using D = Domain<Mesh, Element, NodeIndex>;
    py::class_<D, std::unique_ptr<D, py::nodelete>> c(m, name.c_str());

    c.def("canonical_dimension", &D::canonical_dimension);
    c.def("number_of_nodes_per_elements", &D::number_of_nodes_per_elements);
    c.def("number_of_elements", &D::number_of_elements);
    c.def("element", [](const D & domain, const UNSIGNED_INTEGER_TYPE & element_id) {
        return domain.element(element_id);
    }, py::arg("element_id"));

    c.def("element_indices", &D::element_indices, py::arg("index"));
    c.def("mesh", &D::mesh);

    // Mesh's add_domain binding for Domain<Element, NodeIndex> type
    m.def("add_domain", [](Mesh & mesh, const std::string & domain_name, const Element &, const Eigen::Matrix<NodeIndex, Eigen::Dynamic, geometry::traits<Element>::NumberOfNodesAtCompileTime> & node_indices) {
        return mesh.template add_domain<Element, NodeIndex>(domain_name, node_indices);
    }, py::arg("domain_name"), py::arg("element_type"), py::arg("node_indices"));
    m.def("add_domain", [](Mesh & mesh, const Element &, const Eigen::Matrix<NodeIndex, Eigen::Dynamic, geometry::traits<Element>::NumberOfNodesAtCompileTime> & node_indices) {
        return mesh.template add_domain<Element, NodeIndex>(node_indices);
    }, py::arg("element_type"), py::arg("node_indices"));

    m.def("add_domain", [](Mesh & mesh, const std::string & domain_name, const Element &, const Eigen::Matrix<NodeIndex, Eigen::Dynamic, Eigen::Dynamic> & node_indices) {
        if (geometry::traits<Element>::NumberOfNodesAtCompileTime != caribou::Dynamic and node_indices.cols() != geometry::traits<Element>::NumberOfNodesAtCompileTime) {
            std::ostringstream ss;
            ss << "You tried to create a domain containing elements of " << geometry::traits<Element>::NumberOfNodesAtCompileTime
               << " nodes, but elements having " << node_indices.cols() << " nodes were found instead. The node indices "
               << "matrix should have NxM indices where N is the number of elements, and M is the number of nodes per element.";
            throw py::key_error(ss.str());
        }

        return mesh.template add_domain<Element, NodeIndex>(domain_name, node_indices);
    }, py::arg("domain_name"), py::arg("element_type"), py::arg("node_indices"));

    m.def("add_domain", [](Mesh & mesh, const Element &, const Eigen::Matrix<NodeIndex, Eigen::Dynamic, Eigen::Dynamic> & node_indices) {
        if (geometry::traits<Element>::NumberOfNodesAtCompileTime != caribou::Dynamic and node_indices.cols() != geometry::traits<Element>::NumberOfNodesAtCompileTime) {
            std::ostringstream ss;
            ss << "You tried to create a domain containing elements of " << geometry::traits<Element>::NumberOfNodesAtCompileTime
               << " nodes, but elements having " << node_indices.cols() << " nodes were found instead. The node indices "
               << "matrix should have NxM indices where N is the number of elements, and M is the number of nodes per element.";
            throw py::key_error(ss.str());
        }

        return mesh.template add_domain<Element, NodeIndex>(node_indices);
    }, py::arg("element_type"), py::arg("node_indices"));

    declare_barycentric_container(c);
}

template <typename Mesh, typename Element>
void declare_domain(py::class_<Mesh> & m, const std::string & name) {
//    todo(jnbrunet2000@gmail.com): Uncomment once we allow a no_copy options for the matrices (pass by ref of a np.array)
//    declare_domain<Mesh, Element, int>(m, name+"int");
//    declare_domain<Mesh, Element, long int>(m, name+"lint");
//    declare_domain<Mesh, Element, long long int>(m, name+"llint");
    declare_domain<Mesh, Element, unsigned int>(m, name+"uint");
//    declare_domain<Mesh, Element, unsigned long int>(m, name+"ulint");
//    declare_domain<Mesh, Element, unsigned long long int>(m, name+"ullint");
}

template<typename Mesh>
void declare_domains(py::class_<Mesh> & m) {
    using namespace geometry;

    // Segments
    if constexpr (Mesh::Dimension >= 1) {
        declare_domain<Mesh, Segment<Mesh::Dimension, Linear>>(m, "Domain" + std::to_string(Mesh::Dimension) + "DSegmentLinear");
        declare_domain<Mesh, Segment<Mesh::Dimension, Quadratic>>(m, "Domain" + std::to_string(Mesh::Dimension) + "DSegmentQuadratic");
    }

    if constexpr (Mesh::Dimension >= 2) {
        // Triangles
        declare_domain<Mesh, Triangle<Mesh::Dimension, Linear>>(m, "Domain"+std::to_string(Mesh::Dimension)+"DTriangleLinear");
        declare_domain<Mesh, Triangle<Mesh::Dimension, Quadratic>>(m, "Domain"+std::to_string(Mesh::Dimension)+"DTriangleQuadratic");

        // Quads
        declare_domain<Mesh, Quad<Mesh::Dimension, Linear>>(m, "Domain"+std::to_string(Mesh::Dimension)+"DQuadLinear");
        declare_domain<Mesh, Quad<Mesh::Dimension, Quadratic>>(m, "Domain"+std::to_string(Mesh::Dimension)+"DQuadQuadratic");
    }

    if constexpr (Mesh::Dimension == 3) {
        // Hexahedrons
        declare_domain<Mesh, Hexahedron<Linear>>(m, "Domain3DHexahedronLinear");
        declare_domain<Mesh, Hexahedron<Quadratic>>(m, "Domain3DHexahedronQuadratic");

        // Tetrahedrons
        declare_domain<Mesh, Tetrahedron<Linear>>(m, "Domain3DTetrahedronLinear");
        declare_domain<Mesh, Tetrahedron<Quadratic>>(m, "Domain3DTetrahedronQuadratic");
    }
}
} // namespace caribou::topology::bindings