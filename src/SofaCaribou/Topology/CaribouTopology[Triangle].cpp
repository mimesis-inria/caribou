#include <SofaCaribou/config.h>
#include <SofaCaribou/Topology/CaribouTopology.inl>
#include <SofaCaribou/Topology/CaribouTopology[Triangle].h>

DISABLE_ALL_WARNINGS_BEGIN
#include <sofa/core/ObjectFactory.h>
#include <SofaBaseTopology/TriangleSetTopologyContainer.h>
DISABLE_ALL_WARNINGS_END

using sofa::core::RegisterObject;
using namespace caribou::geometry;
using namespace sofa::component::topology;
using namespace sofa::core::topology;
using namespace caribou;

namespace SofaCaribou::topology {

// Triangle 2D linear specialization
template<>
auto CaribouTopology<Triangle<_2D>>::templateName(const CaribouTopology<Triangle<_2D>> *) -> std::string {
    return "Triangle_2D";
}

template <>
auto CaribouTopology<Triangle<_2D>>::mesh_is_compatible(const BaseMeshTopology * topology) -> bool {
    return (
        dynamic_cast<const TriangleSetTopologyContainer*>(topology) != nullptr
    );
}

template <>
auto CaribouTopology<Triangle<_2D>>::get_indices_data_from(const BaseMeshTopology * topology) -> sofa::core::objectmodel::BaseData * {
    return topology->findData("triangles");
}

// Triangle 3D linear specialization
template<>
auto
CaribouTopology<Triangle<_3D>>::templateName(const CaribouTopology<Triangle<_3D>> *) -> std::string {
    return "Triangle";
}

template <>
auto CaribouTopology<Triangle<_3D>>::mesh_is_compatible(const BaseMeshTopology * topology) -> bool {
    return (
            dynamic_cast<const TriangleSetTopologyContainer*>(topology) != nullptr
    );
}

template <>
auto CaribouTopology<Triangle<_3D>>::get_indices_data_from(const BaseMeshTopology * topology) -> sofa::core::objectmodel::BaseData * {
    return topology->findData("triangles");
}

// This will force the compiler to compile the class with some template type
template class CaribouTopology<Triangle<_2D>>;
template class CaribouTopology<Triangle<_3D>>;

} // namespace SofaCaribou::topology

namespace sofa::core::objectmodel {
using namespace SofaCaribou::topology;

[[maybe_unused]]
static int _c_ = RegisterObject("Caribou topology")
.add<CaribouTopology<Triangle<_2D>>>()
.add<CaribouTopology<Triangle<_3D>>>()
;
}