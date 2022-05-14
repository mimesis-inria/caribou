#include <SofaCaribou/config.h>
#include <SofaCaribou/Topology/CaribouTopology.inl>
#include <SofaCaribou/Topology/CaribouTopology[Hexahedron20].h>

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

// Hexahedron Quadratic specialization
template<>
auto
CaribouTopology<Hexahedron20>::templateName(const CaribouTopology<Hexahedron20> *) -> std::string {
    return "Hexahedron20";
}

template <>
auto CaribouTopology<Hexahedron20>::mesh_is_compatible(const BaseMeshTopology *) -> bool {
    // We cannot create a quadratic topology from a linear one
    return false;
}

template class CaribouTopology<Hexahedron20>;

} // namespace SofaCaribou::topology

namespace sofa::core::objectmodel {
using namespace SofaCaribou::topology;

[[maybe_unused]]
static int _c_ = RegisterObject("Caribou topology")
        .add<CaribouTopology<Hexahedron20>>();
}