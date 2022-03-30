#pragma once

#include <Caribou/Geometry/Hexahedron20.h>

#include <SofaCaribou/Mapping/CaribouBarycentricMapping.h>
#include <SofaCaribou/Topology/CaribouTopology[Hexahedron20].h>

DISABLE_ALL_WARNINGS_BEGIN
#include <sofa/version.h>
#include <sofa/defaulttype/VecTypes.h>
DISABLE_ALL_WARNINGS_END

namespace SofaCaribou::mapping {

// Hexahedron quadratic specialization
extern template class CaribouBarycentricMapping<caribou::geometry::Hexahedron20, sofa::defaulttype::Vec3Types>;

} // namespace SofaCaribou::mapping
