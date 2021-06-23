#include <SofaCaribou/config.h>
#include <SofaCaribou/Topology/CaribouTopology.inl>
#include <SofaCaribou/Topology/CaribouTopology[Hexahedron].h>

DISABLE_ALL_WARNINGS_BEGIN
#include <sofa/core/ObjectFactory.h>
#include <SofaBaseTopology/HexahedronSetTopologyContainer.h>
DISABLE_ALL_WARNINGS_END

using sofa::core::RegisterObject;
using namespace caribou::geometry;
using namespace sofa::component::topology;
using namespace sofa::core::objectmodel;
using namespace sofa::core::topology;
using namespace caribou;

namespace SofaCaribou::topology {

// Hexahedron linear specialization
template<>
auto CaribouTopology<Hexahedron<Linear>>::templateName(const CaribouTopology<Hexahedron<Linear>> *) -> std::string {
    return "Hexahedron";
}

template <>
auto CaribouTopology<Hexahedron<Linear>>::mesh_is_compatible(const BaseMeshTopology * topology) -> bool {
    return (
            dynamic_cast<const HexahedronSetTopologyContainer*>(topology) != nullptr
    );
}

template <>
auto CaribouTopology<Hexahedron<Linear>>::get_indices_data_from(
        const BaseMeshTopology * topology) -> sofa::core::objectmodel::BaseData * {
    return topology->findData("hexahedra");
}

// Hexahedron Quadratic specialization
template<>
auto
CaribouTopology<Hexahedron<Quadratic>>::templateName(const CaribouTopology<Hexahedron<Quadratic>> *) -> std::string {
    return "Hexahedron20";
}

template <>
auto CaribouTopology<Hexahedron<Quadratic>>::mesh_is_compatible(const BaseMeshTopology *) -> bool {
    // We cannot create a quadratic topology from a linear one
    return false;
}

// This will force the compiler to compile the class with some template type
template
class CaribouTopology<Hexahedron<Linear>>;

template
class CaribouTopology<Hexahedron<Quadratic>>;

} // namespace SofaCaribou::topology

namespace sofa::core::objectmodel {
using namespace SofaCaribou::topology;

[[maybe_unused]]
static int _c_ = RegisterObject("Caribou topology")
        .add<CaribouTopology<Hexahedron<Linear>>>()
        .add<CaribouTopology<Hexahedron<Quadratic>>>();
}