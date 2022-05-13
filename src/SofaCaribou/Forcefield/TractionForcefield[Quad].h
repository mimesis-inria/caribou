#pragma once
#include <SofaCaribou/Forcefield/TractionForcefield.h>
#include <SofaCaribou/Forcefield/CaribouForcefield[Quad].h>
#include <Caribou/Geometry/Quad.h>

namespace SofaCaribou::forcefield {

// Quad linear specialization
extern template class TractionForcefield<caribou::geometry::Quad<caribou::_2D>>;
extern template class TractionForcefield<caribou::geometry::Quad<caribou::_3D>>;

}