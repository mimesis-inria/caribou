#pragma once
#include <SofaCaribou/FEniCS/Forcefield/HyperelasticForcefield_FEniCS.h>
#include <SofaCaribou/Forcefield/CaribouForcefield[Hexahedron_FEniCS].h>
#include <Caribou/Geometry/Hexahedron_FEniCS.h>

namespace SofaCaribou::forcefield {

// Hexahedron_FEniCS linear specialization
extern template class HyperelasticForcefield_FEniCS<caribou::geometry::Hexahedron_FEniCS>;

}
