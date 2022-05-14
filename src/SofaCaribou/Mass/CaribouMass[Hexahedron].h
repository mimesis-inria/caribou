#pragma once

#include <SofaCaribou/Mass/CaribouMass.h>
#include <SofaCaribou/Topology/CaribouTopology[Hexahedron].h>
#include <Caribou/Geometry/Hexahedron.h>

namespace SofaCaribou::mass {

// Hexahedron linear specialization
extern template class CaribouMass<caribou::geometry::Hexahedron>;

} // namespace SofaCaribou::mass