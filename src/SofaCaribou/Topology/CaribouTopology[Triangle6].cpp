#include <SofaCaribou/config.h>
#include <SofaCaribou/Topology/CaribouTopology.inl>
#include <SofaCaribou/Topology/CaribouTopology[Triangle6].h>

DISABLE_ALL_WARNINGS_BEGIN
#include <sofa/core/ObjectFactory.h>
#include <sofa/component/topology/container/dynamic/TriangleSetTopologyContainer.h>
DISABLE_ALL_WARNINGS_END

using sofa::core::RegisterObject;
using namespace caribou::geometry;
using namespace sofa::component::topology;
using namespace sofa::core::topology;
using namespace caribou;

namespace SofaCaribou::topology {

// Triangle 2D quadratic specialization
template<>
auto CaribouTopology<Triangle6<_2D>>::GetCustomTemplateName() -> std::string {
    return "Triangle6_2D";
}

template <>
auto CaribouTopology<Triangle6<_2D>>::mesh_is_compatible(const BaseMeshTopology *) -> bool {
    // We cannot create a quadratic topology from a linear one
    return false;
}

// Triangle 3D quadratic specialization
template<>
auto
CaribouTopology<Triangle6<_3D>>::GetCustomTemplateName() -> std::string {
    return "Triangle6";
}

template <>
auto CaribouTopology<Triangle6<_3D>>::mesh_is_compatible(const BaseMeshTopology *) -> bool {
    // We cannot create a quadratic topology from a linear one
    return false;
}

// This will force the compiler to compile the class with some template type
template class CaribouTopology<Triangle6<_2D>>;
template class CaribouTopology<Triangle6<_3D>>;

} // namespace SofaCaribou::topology

namespace sofa::core::objectmodel {
using namespace SofaCaribou::topology;

[[maybe_unused]]
static int _c_ = RegisterObject("Caribou topology")
.add<CaribouTopology<Triangle6<_2D>>>()
.add<CaribouTopology<Triangle6<_3D>>>()
;
}