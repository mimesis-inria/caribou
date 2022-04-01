#include <SofaCaribou/Python/Topology/CaribouTopology.h>

#include <SofaCaribou/Topology/CaribouTopology[Triangle].h>
#include <SofaCaribou/Topology/CaribouTopology[Triangle6].h>
#include <SofaCaribou/Topology/CaribouTopology[Tetrahedron].h>
#include <SofaCaribou/Topology/CaribouTopology[Tetrahedron10].h>
#include <SofaCaribou/Topology/CaribouTopology[Hexahedron].h>
#include <SofaCaribou/Topology/CaribouTopology[Hexahedron20].h>
#include <SofaCaribou/Topology/CaribouTopology[Quad].h>
#include <SofaCaribou/Topology/CaribouTopology[Quad8].h>

#include <pybind11/pybind11.h>

namespace SofaCaribou::topology::python {

void addCaribouTopology(pybind11::module &m) {
    using namespace caribou;
    using namespace caribou::geometry;
    bind_caribou_topology<Triangle<_2D>>(m);
    bind_caribou_topology<Triangle6<_2D>>(m);
    bind_caribou_topology<Triangle<_3D>>(m);
    bind_caribou_topology<Triangle6<_3D>>(m);

    bind_caribou_topology<Quad<_2D>>(m);
    bind_caribou_topology<Quad8<_2D>>(m);
    bind_caribou_topology<Quad<_3D>>(m);
    bind_caribou_topology<Quad8<_3D>>(m);

    bind_caribou_topology<Tetrahedron>(m);
    bind_caribou_topology<Tetrahedron10>(m);

    bind_caribou_topology<Hexahedron>(m);
    bind_caribou_topology<Hexahedron20>(m);
}

} // namespace SofaCaribou::topology::python