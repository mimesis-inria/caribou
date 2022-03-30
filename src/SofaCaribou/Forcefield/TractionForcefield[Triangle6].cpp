#include <SofaCaribou/Forcefield/TractionForcefield[Triangle6].h>
#include <SofaCaribou/Forcefield/TractionForcefield.inl>

DISABLE_ALL_WARNINGS_BEGIN
#include <sofa/core/ObjectFactory.h>
DISABLE_ALL_WARNINGS_END

using sofa::core::RegisterObject;
using namespace caribou::geometry;
using namespace caribou;

namespace SofaCaribou::forcefield {

// -----------------------------------
// Quad quadratic specialization
// -----------------------------------

// This will force the compiler to compile the following templated class
template class TractionForcefield<Triangle6<_2D>>;
template class TractionForcefield<Triangle6<_3D>>;


} // namespace SofaCaribou::forcefield

namespace sofa::core::objectmodel {
using namespace SofaCaribou::forcefield;

[[maybe_unused]]
static int _c_ = RegisterObject("Caribou traction force field")
        .add<TractionForcefield<Triangle6<_2D>>>()
        .add<TractionForcefield<Triangle6<_3D>>>();
}
