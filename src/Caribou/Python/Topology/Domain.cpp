#include <pybind11/pybind11.h>

#include <Caribou/Geometry/Hexahedron.h>
#include <Caribou/Geometry/Hexahedron20.h>
#include <Caribou/Geometry/Quad.h>
#include <Caribou/Geometry/Quad8.h>
#include <Caribou/Geometry/Segment.h>
#include <Caribou/Geometry/Segment3.h>
#include <Caribou/Geometry/Triangle.h>
#include <Caribou/Geometry/Triangle6.h>
#include <Caribou/Geometry/Tetrahedron.h>
#include <Caribou/Geometry/Tetrahedron10.h>
#include <Caribou/Topology/Domain.h>

#include <Caribou/Python/Topology/Domain.h>

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
        declare_domain<Segment<Dimension>>(m, "Domain" + std::to_string(Dimension) + "DSegment");
        declare_domain<Segment3<Dimension>>(m, "Domain" + std::to_string(Dimension) + "DSegment3");
    }

    if constexpr (Dimension >= 2) {
        // Triangles
        declare_domain<Triangle<Dimension>>(m, "Domain" + std::to_string(Dimension) + "DTriangle");
        declare_domain<Triangle6<Dimension>>(m, "Domain" + std::to_string(Dimension) + "DTriangle6");

        // Quads
        declare_domain<Quad<Dimension>>(m, "Domain" + std::to_string(Dimension) + "DQuad");
        declare_domain<Quad8<Dimension>>(m, "Domain" + std::to_string(Dimension) + "DQuad8");
    }

    if constexpr (Dimension == 3) {
        // Hexahedrons
        declare_domain<Hexahedron>(m, "Domain3DHexahedron");
        declare_domain<Hexahedron20>(m, "Domain3DHexahedron20");

        // Tetrahedrons
        declare_domain<Tetrahedron>(m, "Domain3DTetrahedron");
        declare_domain<Tetrahedron10>(m, "Domain3DTetrahedron10");
    }
}

void create_domain(py::module &m) {
    declare_domain<1>(m); // 1D
    declare_domain<2>(m); // 2D
    declare_domain<3>(m); // 3D
}

} // namespace caribou::topology::bindings
