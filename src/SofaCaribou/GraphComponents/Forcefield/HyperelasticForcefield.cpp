#include <SofaCaribou/GraphComponents/Forcefield/HyperelasticForcefield.inl>

#include <Caribou/Geometry/Triangle.h>
#include <Caribou/Geometry/Quad.h>
#include <Caribou/Geometry/Tetrahedron.h>
#include <Caribou/Geometry/Hexahedron.h>

#include <sofa/core/ObjectFactory.h>
#include <SofaBaseTopology/TetrahedronSetTopologyContainer.h>
#include <SofaBaseTopology/HexahedronSetTopologyContainer.h>


using sofa::core::RegisterObject;
using namespace sofa::component::topology;
using namespace caribou::geometry;
namespace SofaCaribou::GraphComponents::forcefield {

// ------------
// Triangle 2D
// ------------
template <>
auto HyperelasticForcefield<Triangle<2, interpolation::Triangle3>>::number_of_elements() const -> std::size_t {
    if (not d_topology_container.empty()) {
        return d_topology_container->getNbTriangles();
    }
    return 0;
}

template <>
auto HyperelasticForcefield<Triangle<2, interpolation::Triangle3>>::mesh_is_compatible(const sofa::core::topology::BaseMeshTopology * topology) -> bool {
    return (
        dynamic_cast<const TriangleSetTopologyContainer*>   (topology) != nullptr and
        dynamic_cast<const TetrahedronSetTopologyContainer*>(topology) == nullptr
    );
}

template <>
auto HyperelasticForcefield<Triangle<2, interpolation::Triangle3>>::get_element_nodes_indices(const std::size_t & element_id) const -> const Index * {
    return &(d_topology_container->getTriangles()[element_id][0]);
}

template <>
auto HyperelasticForcefield<Triangle<2, interpolation::Triangle3>>::templateName(const HyperelasticForcefield<Triangle<2, interpolation::Triangle3>> *) -> std::string {
    return "Triangle";
}

// ------------
// Quad 2D
// ------------
template <>
auto HyperelasticForcefield<Quad<2, interpolation::Quad4>>::number_of_elements() const -> std::size_t {
    if (not d_topology_container.empty()) {
        return d_topology_container->getNbQuads();
    }
    return 0;
}

template <>
auto HyperelasticForcefield<Quad<2, interpolation::Quad4>>::mesh_is_compatible(const sofa::core::topology::BaseMeshTopology * topology) -> bool {
    return (
        dynamic_cast<const QuadSetTopologyContainer*>      (topology) != nullptr and
        dynamic_cast<const HexahedronSetTopologyContainer*>(topology) == nullptr
    );
}

template <>
auto HyperelasticForcefield<Quad<2, interpolation::Quad4>>::get_element_nodes_indices(const std::size_t & element_id) const -> const Index * {
    return &(d_topology_container->getQuads()[element_id][0]);
}

template <>
auto HyperelasticForcefield<Quad<2, interpolation::Quad4>>::templateName(const HyperelasticForcefield<Quad<2, interpolation::Quad4>> *) -> std::string {
    return "Quad";
}

// ------------
// Tetrahedron
// ------------
template <>
auto HyperelasticForcefield<Tetrahedron<interpolation::Tetrahedron4>>::number_of_elements() const -> std::size_t {
    if (not d_topology_container.empty()) {
        return d_topology_container->getNbTetras();
    }
    return 0;
}

template <>
auto HyperelasticForcefield<Tetrahedron<interpolation::Tetrahedron4>>::mesh_is_compatible(const sofa::core::topology::BaseMeshTopology * topology) -> bool {
    return (
        dynamic_cast<const TetrahedronSetTopologyContainer*>(topology) != nullptr
    );
}

template <>
auto HyperelasticForcefield<Tetrahedron<interpolation::Tetrahedron4>>::get_element_nodes_indices(const std::size_t & element_id) const -> const Index * {
    return &(d_topology_container->getTetras()[element_id][0]);
}

template <>
auto HyperelasticForcefield<Tetrahedron<interpolation::Tetrahedron4>>::templateName(const HyperelasticForcefield<Tetrahedron<interpolation::Tetrahedron4>> *) -> std::string {
    return "Tetrahedron";
}

// ------------
// Hexahedron
// ------------
template <>
auto HyperelasticForcefield<Hexahedron<interpolation::Hexahedron8>>::number_of_elements() const -> std::size_t {
    if (not d_topology_container.empty()) {
        return d_topology_container->getNbHexas();
    }
    return 0;
}

template <>
auto HyperelasticForcefield<Hexahedron<interpolation::Hexahedron8>>::mesh_is_compatible(const sofa::core::topology::BaseMeshTopology * topology) -> bool {
    return (
        dynamic_cast<const HexahedronSetTopologyContainer*>(topology) != nullptr
    );
}

template <>
auto HyperelasticForcefield<Hexahedron<interpolation::Hexahedron8>>::get_element_nodes_indices(const std::size_t & element_id) const -> const Index * {
    return &(d_topology_container->getHexas()[element_id][0]);
}

template <>
auto HyperelasticForcefield<Hexahedron<interpolation::Hexahedron8>>::templateName(const HyperelasticForcefield<Hexahedron<interpolation::Hexahedron8>> *) -> std::string {
    return "Hexahedron";
}
}

namespace sofa::core::objectmodel {
using namespace SofaCaribou::GraphComponents::forcefield;

static int HyperelasticForcefieldClass = RegisterObject("Caribou hyperelastic FEM Forcefield")
    .add< HyperelasticForcefield<Tetrahedron<interpolation::Tetrahedron4>> >()
    .add< HyperelasticForcefield<Hexahedron <interpolation::Hexahedron8>> >();
}