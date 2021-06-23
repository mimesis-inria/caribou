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
auto CaribouTopology<Tetrahedron<Linear>>::templateName(const CaribouTopology<Tetrahedron<Linear>> *) -> std::string {
    return "Tetrahedron";
}

template <>
auto CaribouTopology<Tetrahedron<Linear>>::mesh_is_compatible(const BaseMeshTopology * topology) -> bool {
    return (
        dynamic_cast<const TetrahedronSetTopologyContainer*>(topology) != nullptr
    );
}

template <>
auto CaribouTopology<Tetrahedron<Linear>>::get_indices_data_from(
        const BaseMeshTopology * topology) -> sofa::core::objectmodel::BaseData * {
    return topology->findData("tetrahedra");
}

// Tetrahedron Quadratic specialization
template<>
auto
CaribouTopology<Tetrahedron<Quadratic>>::templateName(const CaribouTopology<Tetrahedron<Quadratic>> *) -> std::string {
    return "Tetrahedron10";
}

template <>
auto CaribouTopology<Tetrahedron<Quadratic>>::mesh_is_compatible(const BaseMeshTopology *) -> bool {
    // We cannot create a quadratic topology from a linear one
    return false;
}

// This will force the compiler to compile the class with some template type
template
class CaribouTopology<Tetrahedron<Linear>>;

template
class CaribouTopology<Tetrahedron<Quadratic>>;

} // namespace SofaCaribou::topology

namespace sofa::core::objectmodel {
using namespace SofaCaribou::topology;

[[maybe_unused]]
static int _c_ = RegisterObject("Caribou topology")
        .add<CaribouTopology<Tetrahedron<Linear>>>()
        .add<CaribouTopology<Tetrahedron<Quadratic>>>();
}