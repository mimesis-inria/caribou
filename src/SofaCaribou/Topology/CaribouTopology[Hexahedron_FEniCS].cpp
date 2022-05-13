#include <SofaCaribou/config.h>
#include <SofaCaribou/Topology/CaribouTopology.inl>
#include <SofaCaribou/Topology/CaribouTopology[Hexahedron_FEniCS].h>

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

// Hexahedron_FEniCS linear specialization
template<>
auto CaribouTopology<Hexahedron_FEniCS>::templateName(const CaribouTopology<Hexahedron_FEniCS> *) -> std::string {
    return "Hexahedron_FEniCS";
}

template <>
auto CaribouTopology<Hexahedron_FEniCS>::mesh_is_compatible(const BaseMeshTopology * topology) -> bool {
    return (
            dynamic_cast<const HexahedronSetTopologyContainer*>(topology) != nullptr
    );
}

template <>
auto CaribouTopology<Hexahedron_FEniCS>::get_indices_data_from(
        const BaseMeshTopology * topology) -> sofa::core::objectmodel::BaseData * {
    return topology->findData("hexahedra");
}

// This will force the compiler to compile the class with some template type
template
class CaribouTopology<Hexahedron_FEniCS>;


} // namespace SofaCaribou::topology

namespace sofa::core::objectmodel {
using namespace SofaCaribou::topology;

[[maybe_unused]]
static int _c_ = RegisterObject("Caribou topology")
        .add<CaribouTopology<Hexahedron_FEniCS>>();
}
