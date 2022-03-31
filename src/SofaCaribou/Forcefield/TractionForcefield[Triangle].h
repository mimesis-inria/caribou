#pragma once
#include <SofaCaribou/Forcefield/TractionForcefield.h>
#include <SofaCaribou/Forcefield/CaribouForcefield[Triangle].h>
#include <Caribou/Geometry/Triangle.h>

namespace SofaCaribou::forcefield {

// Triangle linear specialization
extern template class TractionForcefield<caribou::geometry::Triangle<caribou::_2D>>;
extern template class TractionForcefield<caribou::geometry::Triangle<caribou::_3D>>;

}