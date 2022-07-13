#pragma once

#include <SofaCaribou/Mass/CaribouMass.h>
#include <SofaCaribou/Topology/CaribouTopology[Hexahedron_FEniCS20].h>
#include <Caribou/Geometry/Hexahedron_FEniCS20.h>

namespace SofaCaribou::mass {

// Hexahedron quadratic specialization
extern template class CaribouMass<caribou::geometry::Hexahedron_FEniCS20>;

} // namespace SofaCaribou::mass
