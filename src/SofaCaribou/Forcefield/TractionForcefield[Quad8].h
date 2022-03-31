#pragma once
#include <SofaCaribou/Forcefield/TractionForcefield.h>
#include <SofaCaribou/Forcefield/CaribouForcefield[Quad8].h>
#include <Caribou/Geometry/Quad8.h>

namespace SofaCaribou::forcefield {
// Quad quadratic specialization

extern template class TractionForcefield<caribou::geometry::Quad8<caribou::_2D>>;
extern template class TractionForcefield<caribou::geometry::Quad8<caribou::_3D>>;
}