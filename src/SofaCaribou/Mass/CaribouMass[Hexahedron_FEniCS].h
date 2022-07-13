#pragma once

#include <SofaCaribou/Mass/CaribouMass.h>
#include <SofaCaribou/Topology/CaribouTopology[Hexahedron_FEniCS].h>
#include <Caribou/Geometry/Hexahedron_FEniCS.h>

namespace SofaCaribou::mass {

// Hexahedron_FEniCS linear specialization
extern template class CaribouMass<caribou::geometry::Hexahedron_FEniCS>;

} // namespace SofaCaribou::mass
