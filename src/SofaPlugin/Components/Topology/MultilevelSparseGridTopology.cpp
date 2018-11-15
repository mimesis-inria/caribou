#include <sofa/core/ObjectFactory.h>
#include "MultilevelSparseGridTopology.inl"

namespace caribou {
namespace SofaPlugin {
namespace components {
namespace topology {

using sofa::defaulttype::Vec2Types;
using sofa::defaulttype::Vec3Types;



// This will force the compiler to compile the class with some template type
template class MultilevelSparseGridTopology<Vec2Types>;
template class MultilevelSparseGridTopology<Vec3Types>;

// Add the sofa component to the object factory
int MultilevelSparseGridTopologyClass = sofa::core::RegisterObject("Caribou MultilevelSparseGridTopology")
        .add< MultilevelSparseGridTopology<Vec2Types> >()
        .add< MultilevelSparseGridTopology<Vec3Types> >(true)
;

} // namespace topology
} // namespace components
} // namespace SofaPlugin
} // namespace caribou
