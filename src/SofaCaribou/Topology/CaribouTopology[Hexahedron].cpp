#include <SofaCaribou/Topology/CaribouTopology.inl>
#include <SofaCaribou/Topology/CaribouTopology[Hexahedron].h>

#include <sofa/core/ObjectFactory.h>

using sofa::core::RegisterObject;
using namespace caribou::geometry;
using namespace caribou;

namespace SofaCaribou::topology {

// Hexahedron linear specialization
template<>
auto CaribouTopology<Hexahedron<Linear>>::templateName(const CaribouTopology<Hexahedron<Linear>> *) -> std::string {
    return "Hexahedron";
}

// Hexahedron Quadratic specialization
template<>
auto
CaribouTopology<Hexahedron<Quadratic>>::templateName(const CaribouTopology<Hexahedron<Quadratic>> *) -> std::string {
    return "Hexahedron20";
}

// This will force the compiler to compile the class with some template type
template
class CaribouTopology<Hexahedron<Linear>>;

template
class CaribouTopology<Hexahedron<Quadratic>>;

} // namespace SofaCaribou::topology

namespace sofa::core::objectmodel {
using namespace SofaCaribou::topology;

[[maybe_unused]]
static int _c_ = RegisterObject("Caribou topology")
        .add<CaribouTopology<Hexahedron<Linear>>>()
        .add<CaribouTopology<Hexahedron<Quadratic>>>();
}