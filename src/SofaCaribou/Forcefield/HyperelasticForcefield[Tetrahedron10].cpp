#include <SofaCaribou/Forcefield/HyperelasticForcefield[Tetrahedron10].h>
#include <SofaCaribou/Forcefield/HyperelasticForcefield.inl>

DISABLE_ALL_WARNINGS_BEGIN
#include <sofa/core/ObjectFactory.h>
DISABLE_ALL_WARNINGS_END

using sofa::core::RegisterObject;
using namespace caribou::geometry;
using namespace caribou;

namespace SofaCaribou::forcefield {

// -----------------------------------
// Tetrahedron quadratic specialization
// -----------------------------------

// This will force the compiler to compile the following templated class
template class HyperelasticForcefield<Tetrahedron10>;

} // namespace SofaCaribou::forcefield

namespace sofa::core::objectmodel {
using namespace SofaCaribou::forcefield;

[[maybe_unused]]
static int _c_ = RegisterObject("Caribou hyperelastic force field")
    .add<HyperelasticForcefield<Tetrahedron10>>();
}
