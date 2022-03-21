#pragma once
#include <SofaCaribou/ACEgen/Forcefield/HyperelasticForcefield_ACEgen.h>
#include <SofaCaribou/Forcefield/CaribouForcefield[Hexahedron].h>
#include <Caribou/Geometry/Hexahedron.h>

namespace SofaCaribou::forcefield {

// Hexahedron linear specialization
extern template class HyperelasticForcefield_ACEgen<caribou::geometry::Hexahedron < caribou::Linear>>;


}
