#pragma once

#include <SofaCaribou/Mass/CaribouMass.h>
#include <SofaCaribou/Topology/CaribouTopology[Tetrahedron10].h>
#include <Caribou/Geometry/Tetrahedron10.h>

namespace SofaCaribou::mass {

// Tetrahedron quadratic specialization
extern template class CaribouMass<caribou::geometry::Tetrahedron10>;

} // namespace SofaCaribou::mass