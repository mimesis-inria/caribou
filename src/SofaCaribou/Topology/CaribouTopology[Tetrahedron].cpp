#include <SofaCaribou/Topology/CaribouTopology.inl>
#include <SofaCaribou/Topology/CaribouTopology[Tetrahedron].h>

#include <sofa/core/ObjectFactory.h>

using sofa::core::RegisterObject;
using namespace caribou::geometry;
using namespace caribou;

namespace SofaCaribou::topology {

// Tetrahedron linear specialization
template<>
auto CaribouTopology<Tetrahedron<Linear>>::templateName(const CaribouTopology<Tetrahedron<Linear>> *) -> std::string {
    return "Tetrahedron";
}

// Tetrahedron Quadratic specialization
template<>
auto
CaribouTopology<Tetrahedron<Quadratic>>::templateName(const CaribouTopology<Tetrahedron<Quadratic>> *) -> std::string {
    return "Tetrahedron10";
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