#include <SofaCaribou/Python/Topology/CaribouTopology.h>

#include <SofaCaribou/Topology/CaribouTopology[Triangle].h>
#include <SofaCaribou/Topology/CaribouTopology[Tetrahedron].h>
#include <SofaCaribou/Topology/CaribouTopology[Hexahedron].h>
#include <SofaCaribou/Topology/CaribouTopology[Quad].h>

#include <pybind11/pybind11.h>

namespace SofaCaribou::topology::python {

void addCaribouTopology(pybind11::module &m) {
    using namespace caribou;
    using namespace caribou::geometry;
    bind_caribou_topology<Triangle<_2D, Linear>>(m);
    bind_caribou_topology<Triangle<_2D, Quadratic>>(m);
    bind_caribou_topology<Triangle<_3D, Linear>>(m);
    bind_caribou_topology<Triangle<_3D, Quadratic>>(m);

    bind_caribou_topology<Quad<_2D, Linear>>(m);
    bind_caribou_topology<Quad<_2D, Quadratic>>(m);
    bind_caribou_topology<Quad<_3D, Linear>>(m);
    bind_caribou_topology<Quad<_3D, Quadratic>>(m);

    bind_caribou_topology<Tetrahedron<Linear>>(m);
    bind_caribou_topology<Tetrahedron<Quadratic>>(m);

    bind_caribou_topology<Hexahedron<Linear>>(m);
    bind_caribou_topology<Hexahedron<Quadratic>>(m);
}

} // namespace SofaCaribou::topology::python