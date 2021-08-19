#include <SofaCaribou/config.h>
#include <SofaCaribou/Topology/CaribouTopology.inl>
#include <SofaCaribou/Topology/CaribouTopology[Quad].h>

DISABLE_ALL_WARNINGS_BEGIN
#include <sofa/core/ObjectFactory.h>
#include <SofaBaseTopology/QuadSetTopologyContainer.h>
DISABLE_ALL_WARNINGS_END

using sofa::core::RegisterObject;
using namespace caribou::geometry;
using namespace sofa::component::topology;
using namespace sofa::core::topology;
using namespace caribou;

namespace SofaCaribou::topology {

// Quad 2D linear specialization
template<>
auto CaribouTopology<Quad<_2D, Linear>>::templateName(const CaribouTopology<Quad<_2D, Linear>> *) -> std::string {
    return "Quad_2D";
}

template <>
auto CaribouTopology<Quad<_2D, Linear>>::mesh_is_compatible(const BaseMeshTopology * topology) -> bool {
    return (
    dynamic_cast<const QuadSetTopologyContainer*>(topology) != nullptr
    );
}

template <>
auto CaribouTopology<Quad<_2D, Linear>>::get_indices_data_from(const BaseMeshTopology * topology) -> sofa::core::objectmodel::BaseData * {
    return topology->findData("quads");
}

// Quad 2D quadratic specialization
template<>
auto CaribouTopology<Quad<_2D, Quadratic>>::templateName(const CaribouTopology<Quad<_2D, Quadratic>> *) -> std::string {
    return "Quad8_2D";
}

template <>
auto CaribouTopology<Quad<_2D, Quadratic>>::mesh_is_compatible(const BaseMeshTopology *) -> bool {
    // We cannot create a quadratic topology from a linear one
    return false;
}

// Quad 3D linear specialization
template<>
auto
CaribouTopology<Quad<_3D, Linear>>::templateName(const CaribouTopology<Quad<_3D, Linear>> *) -> std::string {
    return "Quad";
}

template <>
auto CaribouTopology<Quad<_3D, Linear>>::mesh_is_compatible(const BaseMeshTopology * topology) -> bool {
    return (
            dynamic_cast<const QuadSetTopologyContainer*>(topology) != nullptr
    );
}

template <>
auto CaribouTopology<Quad<_3D, Linear>>::get_indices_data_from(const BaseMeshTopology * topology) -> sofa::core::objectmodel::BaseData * {
    return topology->findData("quads");
}

// Quad 3D quadratic specialization
template<>
auto
CaribouTopology<Quad<_3D, Quadratic>>::templateName(const CaribouTopology<Quad<_3D, Quadratic>> *) -> std::string {
    return "Quad8";
}

template <>
auto CaribouTopology<Quad<_3D, Quadratic>>::mesh_is_compatible(const BaseMeshTopology *) -> bool {
    // We cannot create a quadratic topology from a linear one
    return false;
}

// This will force the compiler to compile the class with some template type
template class CaribouTopology<Quad<_2D, Linear>>;
template class CaribouTopology<Quad<_2D, Quadratic>>;
template class CaribouTopology<Quad<_3D, Linear>>;
template class CaribouTopology<Quad<_3D, Quadratic>>;

} // namespace SofaCaribou::topology

namespace sofa::core::objectmodel {
using namespace SofaCaribou::topology;

[[maybe_unused]]
static int _c_ = RegisterObject("Caribou topology")
.add<CaribouTopology<Quad<_2D, Linear>>>()
.add<CaribouTopology<Quad<_2D, Quadratic>>>()
.add<CaribouTopology<Quad<_3D, Linear>>>()
.add<CaribouTopology<Quad<_3D, Quadratic>>>()
;
}