#include <SofaCaribou/FEniCS/Forcefield/HyperelasticForcefield_FEniCS[Hexahedron_FEniCS].h>
#include <SofaCaribou/FEniCS/Forcefield/HyperelasticForcefield_FEniCS.inl>

DISABLE_ALL_WARNINGS_BEGIN
#include <sofa/core/ObjectFactory.h>
DISABLE_ALL_WARNINGS_END

using sofa::core::RegisterObject;
using namespace caribou::geometry;
using namespace caribou;

namespace SofaCaribou::forcefield {

// --------------------------------
// Hexahedron_FEniCS linear specialization
// --------------------------------

// This will force the compiler to compile the following templated class
template class HyperelasticForcefield_FEniCS<Hexahedron_FEniCS < Linear>>;

// -----------------------------------
// Hexahedron_FEniCS quadratic specialization
// -----------------------------------

// This will force the compiler to compile the following templated class
template class HyperelasticForcefield_FEniCS<Hexahedron_FEniCS < Quadratic>>;

} // namespace SofaCaribou::forcefield

namespace sofa::core::objectmodel {
using namespace SofaCaribou::forcefield;

[[maybe_unused]]
static int _c_ = RegisterObject("Caribou hyperelastic force field")
    .add<HyperelasticForcefield_FEniCS<Hexahedron_FEniCS<Linear>>>()
    .add<HyperelasticForcefield_FEniCS<Hexahedron_FEniCS<Quadratic>>>();
}
