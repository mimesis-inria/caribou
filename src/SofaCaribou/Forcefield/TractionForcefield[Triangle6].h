#pragma once
#include <SofaCaribou/Forcefield/TractionForcefield.h>
#include <SofaCaribou/Forcefield/CaribouForcefield[Triangle6].h>
#include <Caribou/Geometry/Triangle6.h>

namespace SofaCaribou::forcefield {
// Triangle quadratic specialization

extern template class TractionForcefield<caribou::geometry::Triangle6<caribou::_2D>>;
extern template class TractionForcefield<caribou::geometry::Triangle6<caribou::_3D>>;
}