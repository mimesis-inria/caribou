#include <sofa/core/ObjectFactory.h>
#include "FictitiousGrid.inl"

namespace SofaCaribou::GraphComponents::topology {

using sofa::defaulttype::Vec2Types;
using sofa::defaulttype::Vec3Types;

// This will force the compiler to compile the class with some template type
template class FictitiousGrid<Vec2Types>;
template class FictitiousGrid<Vec3Types>;

// Add the sofa component to the object factory
int MultilevelSparseGridTopologyClass = sofa::core::RegisterObject("Caribou FictitiousGrid")
        .add< FictitiousGrid<Vec2Types> >()
        .add< FictitiousGrid<Vec3Types> >(true)
;

} // namespace SofaCaribou::GraphComponents::topology
