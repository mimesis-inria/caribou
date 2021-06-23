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
auto CaribouTopology<Triangle<_2D, Linear>>::templateName(const CaribouTopology<Triangle<_2D, Linear>> *) -> std::string {
    return "Triangle_2D";
}

template <>
auto CaribouTopology<Triangle<_2D, Linear>>::mesh_is_compatible(const BaseMeshTopology * topology) -> bool {
    return (
        dynamic_cast<const TriangleSetTopologyContainer*>(topology) != nullptr
    );
}

template <>
auto CaribouTopology<Triangle<_2D, Linear>>::get_indices_data_from(const BaseMeshTopology * topology) -> sofa::core::objectmodel::BaseData * {
    return topology->findData("triangles");
}

// Triangle 2D quadratic specialization
template<>
auto CaribouTopology<Triangle<_2D, Quadratic>>::templateName(const CaribouTopology<Triangle<_2D, Quadratic>> *) -> std::string {
    return "Triangle6_2D";
}

template <>
auto CaribouTopology<Triangle<_2D, Quadratic>>::mesh_is_compatible(const BaseMeshTopology *) -> bool {
    // We cannot create a quadratic topology from a linear one
    return false;
}

// Triangle 3D linear specialization
template<>
auto
CaribouTopology<Triangle<_3D, Linear>>::templateName(const CaribouTopology<Triangle<_3D, Linear>> *) -> std::string {
    return "Triangle";
}

template <>
auto CaribouTopology<Triangle<_3D, Linear>>::mesh_is_compatible(const BaseMeshTopology * topology) -> bool {
    return (
            dynamic_cast<const TriangleSetTopologyContainer*>(topology) != nullptr
    );
}

template <>
auto CaribouTopology<Triangle<_3D, Linear>>::get_indices_data_from(const BaseMeshTopology * topology) -> sofa::core::objectmodel::BaseData * {
    return topology->findData("triangles");
}

// Triangle 3D quadratic specialization
template<>
auto
CaribouTopology<Triangle<_3D, Quadratic>>::templateName(const CaribouTopology<Triangle<_3D, Quadratic>> *) -> std::string {
    return "Triangle6";
}

template <>
auto CaribouTopology<Triangle<_3D, Quadratic>>::mesh_is_compatible(const BaseMeshTopology *) -> bool {
    // We cannot create a quadratic topology from a linear one
    return false;
}

// This will force the compiler to compile the class with some template type
template class CaribouTopology<Triangle<_2D, Linear>>;
template class CaribouTopology<Triangle<_2D, Quadratic>>;
template class CaribouTopology<Triangle<_3D, Linear>>;
template class CaribouTopology<Triangle<_3D, Quadratic>>;

} // namespace SofaCaribou::topology

namespace sofa::core::objectmodel {
using namespace SofaCaribou::topology;

[[maybe_unused]]
static int _c_ = RegisterObject("Caribou topology")
.add<CaribouTopology<Triangle<_2D, Linear>>>()
.add<CaribouTopology<Triangle<_2D, Quadratic>>>()
.add<CaribouTopology<Triangle<_3D, Linear>>>()
.add<CaribouTopology<Triangle<_3D, Quadratic>>>()
;
}