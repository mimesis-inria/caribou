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
auto CaribouTopology<Hexahedron_FEniCS<Linear>>::templateName(const CaribouTopology<Hexahedron_FEniCS<Linear>> *) -> std::string {
    return "Hexahedron_FEniCS";
}

template <>
auto CaribouTopology<Hexahedron_FEniCS<Linear>>::mesh_is_compatible(const BaseMeshTopology * topology) -> bool {
    return (
            dynamic_cast<const HexahedronSetTopologyContainer*>(topology) != nullptr
    );
}

template <>
auto CaribouTopology<Hexahedron_FEniCS<Linear>>::get_indices_data_from(
        const BaseMeshTopology * topology) -> sofa::core::objectmodel::BaseData * {
    return topology->findData("hexahedra");
}

// Hexahedron_FEniCS Quadratic specialization
template<>
auto
CaribouTopology<Hexahedron_FEniCS<Quadratic>>::templateName(const CaribouTopology<Hexahedron_FEniCS<Quadratic>> *) -> std::string {
    return "Hexahedron_FEniCS20";
}

template <>
auto CaribouTopology<Hexahedron_FEniCS<Quadratic>>::mesh_is_compatible(const BaseMeshTopology *) -> bool {
    // We cannot create a quadratic topology from a linear one
    return false;
}

// This will force the compiler to compile the class with some template type
template
class CaribouTopology<Hexahedron_FEniCS<Linear>>;

template
class CaribouTopology<Hexahedron_FEniCS<Quadratic>>;

} // namespace SofaCaribou::topology

namespace sofa::core::objectmodel {
using namespace SofaCaribou::topology;

[[maybe_unused]]
static int _c_ = RegisterObject("Caribou topology")
        .add<CaribouTopology<Hexahedron_FEniCS<Linear>>>()
        .add<CaribouTopology<Hexahedron_FEniCS<Quadratic>>>();
}
