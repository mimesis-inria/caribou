#pragma once
#include <SofaCaribou/Forcefield/HyperelasticForcefield.h>
#include <SofaCaribou/Forcefield/CaribouForcefield[Tetrahedron10].h>
#include <Caribou/Geometry/Tetrahedron10.h>

namespace SofaCaribou::forcefield {

// Tetrahedron quadratic specialization
extern template class HyperelasticForcefield<caribou::geometry::Tetrahedron10>;

}
