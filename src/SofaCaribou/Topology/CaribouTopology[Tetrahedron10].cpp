#include <SofaCaribou/config.h>
#include <SofaCaribou/Topology/CaribouTopology.inl>
#include <SofaCaribou/Topology/CaribouTopology[Tetrahedron10].h>

DISABLE_ALL_WARNINGS_BEGIN
#include <sofa/core/ObjectFactory.h>
#include <sofa/component/topology/container/dynamic/TetrahedronSetTopologyContainer.h>
DISABLE_ALL_WARNINGS_END

using sofa::core::RegisterObject;
using namespace caribou::geometry;
using namespace sofa::component::topology;
using namespace sofa::core::topology;
using namespace caribou;

namespace SofaCaribou::topology {

// Tetrahedron Quadratic specialization
template<>
auto
CaribouTopology<Tetrahedron10>::GetCustomTemplateName() -> std::string {
    return "Tetrahedron10";
}

template <>
auto CaribouTopology<Tetrahedron10>::mesh_is_compatible(const BaseMeshTopology *) -> bool {
    // We cannot create a quadratic topology from a linear one
    return false;
}

// This will force the compiler to compile the class with some template type
template
class CaribouTopology<Tetrahedron10>;

} // namespace SofaCaribou::topology

namespace sofa::core::objectmodel {
using namespace SofaCaribou::topology;

[[maybe_unused]]
static int _c_ = RegisterObject("Caribou topology")
        .add<CaribouTopology<Tetrahedron10>>();
}