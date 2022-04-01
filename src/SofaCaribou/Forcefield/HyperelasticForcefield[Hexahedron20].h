#pragma once
#include <SofaCaribou/Forcefield/HyperelasticForcefield.h>
#include <SofaCaribou/Forcefield/CaribouForcefield[Hexahedron20].h>
#include <Caribou/Geometry/Hexahedron20.h>

namespace SofaCaribou::forcefield {

// Hexahedron quadratic specialization
extern template class HyperelasticForcefield<caribou::geometry::Hexahedron20>;

}