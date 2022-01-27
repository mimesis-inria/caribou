#pragma once
#include <SofaCaribou//FEniCS/Forcefield/HyperelasticForcefield_FEniCS.h>
#include <SofaCaribou/Forcefield/CaribouForcefield[Hexahedron].h>
#include <Caribou/Geometry/Hexahedron.h>

namespace SofaCaribou::forcefield {

// Hexahedron linear specialization
extern template class HyperelasticForcefield_FEniCS<caribou::geometry::Hexahedron < caribou::Linear>>;

// Hexahedron quadratic specialization
extern template class HyperelasticForcefield_FEniCS<caribou::geometry::Hexahedron < caribou::Quadratic>>;

}
