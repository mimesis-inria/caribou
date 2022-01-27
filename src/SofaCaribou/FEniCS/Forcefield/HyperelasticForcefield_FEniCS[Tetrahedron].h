#pragma once
#include <SofaCaribou/FEniCS/Forcefield/HyperelasticForcefield_FEniCS.h>
#include <SofaCaribou/Forcefield/CaribouForcefield[Tetrahedron].h>
#include <Caribou/Geometry/Tetrahedron.h>

namespace SofaCaribou::forcefield {

// Tetrahedron linear specialization
extern template class HyperelasticForcefield_FEniCS<caribou::geometry::Tetrahedron < caribou::Linear>>;

// Tetrahedron quadratic specialization
extern template class HyperelasticForcefield_FEniCS<caribou::geometry::Tetrahedron < caribou::Quadratic>>;

}
