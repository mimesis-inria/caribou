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

#include <Caribou/Python/Caribou.h>
#include <Caribou/Topology/Python/BarycentricContainer.h>

namespace py = pybind11;

namespace caribou::topology::python {
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

    declare_barycentric_container(c);
}

template <typename Mesh, typename Element>
void declare_domain(py::class_<Mesh> & m, const std::string & name) {
    declare_domain<Mesh, Element, int>(m, name+"int");
    declare_domain<Mesh, Element, long int>(m, name+"lint");
    declare_domain<Mesh, Element, long long int>(m, name+"llint");
    declare_domain<Mesh, Element, unsigned int>(m, name+"uint");
    declare_domain<Mesh, Element, unsigned long int>(m, name+"ulint");
    declare_domain<Mesh, Element, unsigned long long int>(m, name+"ullint");
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
} // namespace caribou::topology::python