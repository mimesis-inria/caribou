#include <SofaCaribou/config.h>
#include <SofaCaribou/Topology/CaribouTopology.inl>
#include <SofaCaribou/Topology/CaribouTopology[Tetrahedron].h>

DISABLE_ALL_WARNINGS_BEGIN
#include <sofa/core/ObjectFactory.h>
#include <SofaBaseTopology/TetrahedronSetTopologyContainer.h>
DISABLE_ALL_WARNINGS_END

using sofa::core::RegisterObject;
using namespace caribou::geometry;
using namespace sofa::component::topology;
using namespace sofa::core::topology;
using namespace caribou;

namespace SofaCaribou::topology {

// Tetrahedron linear specialization
template<>
auto CaribouTopology<Tetrahedron>::templateName(const CaribouTopology<Tetrahedron> *) -> std::string {
    return "Tetrahedron";
}

template <>
auto CaribouTopology<Tetrahedron>::mesh_is_compatible(const BaseMeshTopology * topology) -> bool {
    return (
        dynamic_cast<const TetrahedronSetTopologyContainer*>(topology) != nullptr
    );
}

template <>
auto CaribouTopology<Tetrahedron>::get_indices_data_from(
        const BaseMeshTopology * topology) -> sofa::core::objectmodel::BaseData * {
    return topology->findData("tetrahedra");
}

// This will force the compiler to compile the class with some template type
template
class CaribouTopology<Tetrahedron>;

} // namespace SofaCaribou::topology

namespace sofa::core::objectmodel {
using namespace SofaCaribou::topology;

[[maybe_unused]]
static int _c_ = RegisterObject("Caribou topology")
        .add<CaribouTopology<Tetrahedron>>();
}