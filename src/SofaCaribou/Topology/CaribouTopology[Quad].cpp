#include <SofaCaribou/Topology/CaribouTopology.inl>
#include <SofaCaribou/Topology/CaribouTopology[Quad].h>

#include <sofa/core/ObjectFactory.h>

using sofa::core::RegisterObject;
using namespace caribou::geometry;
using namespace caribou;

namespace SofaCaribou::topology {

// Quad 2D linear specialization
template<>
auto CaribouTopology<Quad<_2D, Linear>>::templateName(const CaribouTopology<Quad<_2D, Linear>> *) -> std::string {
    return "Quad_2D";
}

// Quad 2D quadratic specialization
template<>
auto CaribouTopology<Quad<_2D, Quadratic>>::templateName(const CaribouTopology<Quad<_2D, Quadratic>> *) -> std::string {
    return "Quad8_2D";
}

// Quad 3D linear specialization
template<>
auto
CaribouTopology<Quad<_3D, Linear>>::templateName(const CaribouTopology<Quad<_3D, Linear>> *) -> std::string {
    return "Quad";
}

// Quad 3D quadratic specialization
template<>
auto
CaribouTopology<Quad<_3D, Quadratic>>::templateName(const CaribouTopology<Quad<_3D, Quadratic>> *) -> std::string {
    return "Quad8";
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