#pragma once
#include <SofaCaribou/FEniCS/Forcefield/HyperelasticForcefield_FEniCS.h>
#include <SofaCaribou/Forcefield/CaribouForcefield[Hexahedron_FEniCS20].h>
#include <Caribou/Geometry/Hexahedron_FEniCS20.h>

namespace SofaCaribou::forcefield {

// Hexahedron_FEniCS20 linear specialization
extern template class HyperelasticForcefield_FEniCS<caribou::geometry::Hexahedron_FEniCS20>;

}
