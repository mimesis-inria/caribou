#include <SofaCaribou/Mapping/CaribouBarycentricMapping.inl>
#include <SofaCaribou/Mapping/CaribouBarycentricMapping[Hexahedron].h>

DISABLE_ALL_WARNINGS_BEGIN
#include <sofa/core/ObjectFactory.h>
DISABLE_ALL_WARNINGS_END

using sofa::core::RegisterObject;
using namespace caribou::geometry;
using namespace caribou;

namespace SofaCaribou::mapping {

// ---------------------------------
// Hexahedron linear specialization
// ---------------------------------

// This will force the compiler to compile the following templated class
template class CaribouBarycentricMapping<caribou::geometry::Hexahedron<caribou::Linear>, sofa::defaulttype::Vec3Types>;

// -----------------------------------
// Hexahedron quadratic specialization
// -----------------------------------

// This will force the compiler to compile the following templated class
template class CaribouBarycentricMapping<caribou::geometry::Hexahedron<caribou::Quadratic>, sofa::defaulttype::Vec3Types>;

} // namespace SofaCaribou::mapping

namespace sofa::core::objectmodel {
using namespace SofaCaribou::mapping;

[[maybe_unused]]
static int _c_ = RegisterObject("Caribou barycentric mapping")
.add<CaribouBarycentricMapping<Hexahedron<Linear>, sofa::defaulttype::Vec3Types>>()
.add<CaribouBarycentricMapping<Hexahedron<Quadratic>, sofa::defaulttype::Vec3Types>>();
}