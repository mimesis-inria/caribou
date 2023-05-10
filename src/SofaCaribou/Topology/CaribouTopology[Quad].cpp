#include <SofaCaribou/config.h>
#include <SofaCaribou/Topology/CaribouTopology.inl>
#include <SofaCaribou/Topology/CaribouTopology[Quad].h>

DISABLE_ALL_WARNINGS_BEGIN
#include <sofa/core/ObjectFactory.h>
#include <sofa/component/topology/container/dynamic/QuadSetTopologyContainer.h>
DISABLE_ALL_WARNINGS_END

using sofa::core::RegisterObject;
using namespace caribou::geometry;
using namespace sofa::component::topology;
using namespace sofa::core::topology;
using namespace sofa::component::topology::container::dynamic;
using namespace caribou;

namespace SofaCaribou::topology {

// Quad 2D linear specialization
template<>
auto CaribouTopology<Quad<_2D>>::GetCustomTemplateName() -> std::string {
    return "Quad_2D";
}

template <>
auto CaribouTopology<Quad<_2D>>::mesh_is_compatible(const BaseMeshTopology * topology) -> bool {
    return (
    dynamic_cast<const QuadSetTopologyContainer*>(topology) != nullptr
    );
}

template <>
auto CaribouTopology<Quad<_2D>>::get_indices_data_from(const BaseMeshTopology * topology) -> sofa::core::objectmodel::BaseData * {
    return topology->findData("quads");
}

// Quad 3D linear specialization
template<>
auto
CaribouTopology<Quad<_3D>>::GetCustomTemplateName() -> std::string {
    return "Quad";
}

template <>
auto CaribouTopology<Quad<_3D>>::mesh_is_compatible(const BaseMeshTopology * topology) -> bool {
    return (
            dynamic_cast<const QuadSetTopologyContainer*>(topology) != nullptr
    );
}

template <>
auto CaribouTopology<Quad<_3D>>::get_indices_data_from(const BaseMeshTopology * topology) -> sofa::core::objectmodel::BaseData * {
    return topology->findData("quads");
}

// This will force the compiler to compile the class with some template type
template class CaribouTopology<Quad<_2D>>;
template class CaribouTopology<Quad<_3D>>;

} // namespace SofaCaribou::topology

namespace sofa::core::objectmodel {
using namespace SofaCaribou::topology;

[[maybe_unused]]
static int _c_ = RegisterObject("Caribou topology")
.add<CaribouTopology<Quad<_2D>>>()
.add<CaribouTopology<Quad<_3D>>>()
;
}