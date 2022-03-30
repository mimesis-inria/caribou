#pragma once

#include <Caribou/Geometry/Hexahedron.h>

#include <SofaCaribou/Mapping/CaribouBarycentricMapping.h>
#include <SofaCaribou/Topology/CaribouTopology[Hexahedron].h>

DISABLE_ALL_WARNINGS_BEGIN
#include <sofa/version.h>
#include <sofa/defaulttype/VecTypes.h>
DISABLE_ALL_WARNINGS_END

namespace SofaCaribou::mapping {

// Hexahedron linear specialization
extern template class CaribouBarycentricMapping<caribou::geometry::Hexahedron, sofa::defaulttype::Vec3Types>;

} // namespace SofaCaribou::mapping
