#include <SofaCaribou/ACEgen/Forcefield/HyperelasticForcefield_ACEgen[Hexahedron].h>
#include <SofaCaribou/ACEgen/Forcefield/HyperelasticForcefield_ACEgen.inl>

DISABLE_ALL_WARNINGS_BEGIN
#include <sofa/core/ObjectFactory.h>
DISABLE_ALL_WARNINGS_END

using sofa::core::RegisterObject;
using namespace caribou::geometry;
using namespace caribou;

namespace SofaCaribou::forcefield {

// --------------------------------
// Hexahedron linear specialization
// --------------------------------

// This will force the compiler to compile the following templated class
template class HyperelasticForcefield_ACEgen<Hexahedron < Linear>>;

} // namespace SofaCaribou::forcefield

namespace sofa::core::objectmodel {
using namespace SofaCaribou::forcefield;

[[maybe_unused]]
static int _c_ = RegisterObject("Caribou hyperelastic force field")
    .add<HyperelasticForcefield_ACEgen<Hexahedron<Linear>>>();
}
