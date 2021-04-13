#include <pybind11/pybind11.h>

#include <Caribou/Geometry/Hexahedron.h>
#include <Caribou/Geometry/Quad.h>
#include <Caribou/Geometry/Segment.h>
#include <Caribou/Geometry/Triangle.h>
#include <Caribou/Geometry/Tetrahedron.h>
#include <Caribou/Topology/Domain.h>

#include <Caribou/Bindings/Topology/Domain.h>

namespace py = pybind11;

namespace caribou::topology::bindings {

template<typename Element>
void declare_domain(py::module &m, const std::string &name) {
//    todo(jnbrunet2000@gmail.com): Uncomment once we allow a no_copy options for the matrices (pass by ref of a np.array)
//    declare_domain<Element, int>(m, name+"int");
//    declare_domain<Element, long int>(m, name+"lint");
//    declare_domain<Element, long long int>(m, name+"llint");
    declare_domain<Element, unsigned int>(m, name + "uint");
//    declare_domain<Element, unsigned long int>(m, name+"ulint");
//    declare_domain<Element, unsigned long long int>(m, name+"ullint");
}

template<unsigned int Dimension>
void declare_domain(py::module &m) {
    using namespace caribou::geometry;

    // Segments
    if constexpr (Dimension >= 1) {
        declare_domain<Segment<Dimension, Linear>>(m, "Domain" + std::to_string(Dimension) + "DSegmentLinear");
        declare_domain<Segment<Dimension, Quadratic>>(m, "Domain" + std::to_string(Dimension) + "DSegmentQuadratic");
    }

    if constexpr (Dimension >= 2) {
        // Triangles
        declare_domain<Triangle<Dimension, Linear>>(m, "Domain" + std::to_string(Dimension) + "DTriangleLinear");
        declare_domain<Triangle<Dimension, Quadratic>>(m, "Domain" + std::to_string(Dimension) + "DTriangleQuadratic");

        // Quads
        declare_domain<Quad<Dimension, Linear>>(m, "Domain" + std::to_string(Dimension) + "DQuadLinear");
        declare_domain<Quad<Dimension, Quadratic>>(m, "Domain" + std::to_string(Dimension) + "DQuadQuadratic");
    }

    if constexpr (Dimension == 3) {
        // Hexahedrons
        declare_domain<Hexahedron<Linear>>(m, "Domain3DHexahedronLinear");
        declare_domain<Hexahedron<Quadratic>>(m, "Domain3DHexahedronQuadratic");

        // Tetrahedrons
        declare_domain<Tetrahedron<Linear>>(m, "Domain3DTetrahedronLinear");
        declare_domain<Tetrahedron<Quadratic>>(m, "Domain3DTetrahedronQuadratic");
    }
}

void create_domain(py::module &m) {
    declare_domain<1>(m); // 1D
    declare_domain<2>(m); // 2D
    declare_domain<3>(m); // 3D
}

} // namespace caribou::topology::bindings
