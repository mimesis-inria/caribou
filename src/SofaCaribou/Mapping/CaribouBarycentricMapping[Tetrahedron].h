#pragma once

#include <Caribou/Geometry/Tetrahedron.h>

#include <SofaCaribou/Mapping/CaribouBarycentricMapping.h>
#include <SofaCaribou/Topology/CaribouTopology[Tetrahedron].h>

DISABLE_ALL_WARNINGS_BEGIN
#include <sofa/version.h>
#include <sofa/defaulttype/VecTypes.h>
DISABLE_ALL_WARNINGS_END

namespace SofaCaribou::mapping {

// Tetrahedron linear specialization
extern template class CaribouBarycentricMapping<caribou::geometry::Tetrahedron, sofa::defaulttype::Vec3Types>;

} // namespace SofaCaribou::mapping
