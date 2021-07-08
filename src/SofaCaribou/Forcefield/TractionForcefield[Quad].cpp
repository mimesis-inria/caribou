#include <SofaCaribou/Forcefield/TractionForcefield[Quad].h>
#include <SofaCaribou/Forcefield/TractionForcefield.inl>

DISABLE_ALL_WARNINGS_BEGIN
#include <sofa/core/ObjectFactory.h>
DISABLE_ALL_WARNINGS_END

using sofa::core::RegisterObject;
using namespace caribou::geometry;
using namespace caribou;

namespace SofaCaribou::forcefield {

// --------------------------------
// Quad linear specialization
// --------------------------------

// This will force the compiler to compile the following templated class
template class TractionForcefield<Quad <_2D, Linear>>;
template class TractionForcefield<Quad <_3D, Linear>>;

// -----------------------------------
// Quad quadratic specialization
// -----------------------------------

// This will force the compiler to compile the following templated class
template class TractionForcefield<Quad <_2D, Quadratic>>;
template class TractionForcefield<Quad <_3D, Quadratic>>;


} // namespace SofaCaribou::forcefield

namespace sofa::core::objectmodel {
using namespace SofaCaribou::forcefield;

[[maybe_unused]]
static int _c_ = RegisterObject("Caribou traction force field")
        .add<TractionForcefield<Quad <_2D, Linear>>>()
        .add<TractionForcefield<Quad <_3D, Linear>>>()
        .add<TractionForcefield<Quad <_2D, Quadratic>>>()
        .add<TractionForcefield<Quad <_3D, Quadratic>>>();
}
