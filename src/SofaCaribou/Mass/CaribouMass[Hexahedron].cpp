#include <SofaCaribou/Mass/CaribouMass[Hexahedron].h>
#include <SofaCaribou/Topology/CaribouTopology[Hexahedron].h>
#include <SofaCaribou/Mass/CaribouMass.inl>

DISABLE_ALL_WARNINGS_BEGIN
#include <sofa/core/ObjectFactory.h>
DISABLE_ALL_WARNINGS_END

using sofa::core::RegisterObject;
using namespace sofa::core::objectmodel;
using namespace caribou::geometry;
using namespace caribou;

namespace SofaCaribou::mass {

// This will force the compiler to compile the following templated class
template class CaribouMass<Hexahedron<Linear>>;
template class CaribouMass<Hexahedron<Quadratic>>;

[[maybe_unused]]
static int _c_ = RegisterObject("Caribou mass")
        .add<CaribouMass<Hexahedron<Linear>>>()
        .add<CaribouMass<Hexahedron<Quadratic>>>();

} // namespace SofaCaribou::mass