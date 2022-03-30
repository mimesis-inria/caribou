#pragma once

#include <SofaCaribou/Mass/CaribouMass.h>
#include <SofaCaribou/Topology/CaribouTopology[Hexahedron20].h>
#include <Caribou/Geometry/Hexahedron20.h>

namespace SofaCaribou::mass {

// Hexahedron quadratic specialization
extern template class CaribouMass<caribou::geometry::Hexahedron20>;

} // namespace SofaCaribou::mass