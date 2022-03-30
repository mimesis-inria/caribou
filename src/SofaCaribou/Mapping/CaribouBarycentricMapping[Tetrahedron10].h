#pragma once

#include <Caribou/Geometry/Tetrahedron10.h>

#include <SofaCaribou/Mapping/CaribouBarycentricMapping.h>
#include <SofaCaribou/Topology/CaribouTopology[Tetrahedron10].h>

DISABLE_ALL_WARNINGS_BEGIN
#include <sofa/version.h>
#include <sofa/defaulttype/VecTypes.h>
DISABLE_ALL_WARNINGS_END

namespace SofaCaribou::mapping {

// Tetrahedron quadratic specialization
extern template class CaribouBarycentricMapping<caribou::geometry::Tetrahedron10, sofa::defaulttype::Vec3Types>;

} // namespace SofaCaribou::mapping
