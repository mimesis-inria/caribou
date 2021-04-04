#include <SofaCaribou/Forcefield/HyperelasticForcefield.inl>
#include <sofa/core/ObjectFactory.h>
#include <SofaBaseTopology/TetrahedronSetTopologyContainer.h>
#include <SofaBaseTopology/HexahedronSetTopologyContainer.h>


using sofa::core::RegisterObject;
using namespace sofa::component::topology;
using namespace caribou::geometry;
using namespace caribou;
namespace SofaCaribou::forcefield {

// ------------
// Tetrahedron
// ------------
template <>
void HyperelasticForcefield<Tetrahedron<Linear>>::create_domain_from(sofa::core::topology::BaseMeshTopology * topology) {
    // A mesh is required in order to create the domain
    if (not p_mesh) {
        return;
    }

    // Get a reference to the list of tetrahedrons
    const auto & tetrahedrons = topology->getTetrahedra();
    const auto * data = tetrahedrons.data()->data();

    // Create the Domain instance
    using PointID = sofa::core::topology::Topology::PointID;
    p_domain = p_mesh->add_domain<Tetrahedron<Linear>, PointID> (
        data, tetrahedrons.size(), Tetrahedron<Linear>::NumberOfNodesAtCompileTime
    );
}

template <>
auto HyperelasticForcefield<Tetrahedron<Linear>>::mesh_is_compatible(const sofa::core::topology::BaseMeshTopology * topology) -> bool {
    return (
        dynamic_cast<const TetrahedronSetTopologyContainer*>(topology) != nullptr
    );
}

template <>
auto HyperelasticForcefield<Tetrahedron<Linear>>::templateName(const HyperelasticForcefield<Tetrahedron<Linear>> *) -> std::string {
    return "Tetrahedron";
}

// ------------
// Hexahedron
// ------------
template <>
void HyperelasticForcefield<Hexahedron<Linear>>::create_domain_from(sofa::core::topology::BaseMeshTopology * topology) {
    // A mesh is required in order to create the domain
    if (not p_mesh) {
        return;
    }

    // Get a reference to the list of hexahedrons
    const auto & hexahedrons = topology->getHexahedra();
    const auto * data = hexahedrons.data()->data();

    // Create the Domain instance
    using PointID = sofa::core::topology::Topology::PointID;
    p_domain = p_mesh->add_domain<Hexahedron<Linear>, PointID> (
            data, hexahedrons.size(), Hexahedron<Linear>::NumberOfNodesAtCompileTime
    );
}

template <>
auto HyperelasticForcefield<Hexahedron<Linear>>::mesh_is_compatible(const sofa::core::topology::BaseMeshTopology * topology) -> bool {
    return (
        dynamic_cast<const HexahedronSetTopologyContainer*>(topology) != nullptr
    );
}

template <>
auto HyperelasticForcefield<Hexahedron<Linear>>::templateName(const HyperelasticForcefield<Hexahedron<Linear>> *) -> std::string {
    return "Hexahedron";
}


// This will force the compiler to compile the class with some template type
template class HyperelasticForcefield<caribou::geometry::Tetrahedron < caribou::Linear>>;
template class HyperelasticForcefield<caribou::geometry::Hexahedron  < caribou::Linear>>;

}

namespace sofa::core::objectmodel {
using namespace SofaCaribou::forcefield;

static int HyperelasticForcefieldClass = RegisterObject("Caribou hyperelastic FEM Forcefield")
    .add< HyperelasticForcefield<Tetrahedron<Linear>> >()
    .add< HyperelasticForcefield<Hexahedron <Linear>> >();
}
