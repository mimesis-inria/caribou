#include <SofaCaribou/Mass/CaribouMass[Tetrahedron].h>
#include <SofaCaribou/Topology/CaribouTopology[Tetrahedron].h>
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
template class CaribouMass<Tetrahedron<Linear>>;
template class CaribouMass<Tetrahedron<Quadratic>>;

[[maybe_unused]]
static int _c_ = RegisterObject("Caribou mass")
        .add<CaribouMass<Tetrahedron<Linear>>>()
        .add<CaribouMass<Tetrahedron<Quadratic>>>();

} // namespace SofaCaribou::mass