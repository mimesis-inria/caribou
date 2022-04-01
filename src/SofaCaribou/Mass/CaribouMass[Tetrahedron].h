#pragma once

#include <SofaCaribou/Mass/CaribouMass.h>
#include <SofaCaribou/Topology/CaribouTopology[Tetrahedron].h>
#include <Caribou/Geometry/Tetrahedron.h>

namespace SofaCaribou::mass {

// Tetrahedron linear specialization
extern template class CaribouMass<caribou::geometry::Tetrahedron>;

} // namespace SofaCaribou::mass