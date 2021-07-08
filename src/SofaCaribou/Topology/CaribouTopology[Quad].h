#pragma once

#include <SofaCaribou/Topology/CaribouTopology.h>
#include <Caribou/Geometry/Quad.h>

namespace SofaCaribou::topology {

// Quad 2D linear specialization
template<>
auto CaribouTopology<caribou::geometry::Quad<caribou::_2D, caribou::Linear>>::templateName(
        const CaribouTopology<caribou::geometry::Quad<caribou::_2D, caribou::Linear>> *) -> std::string;

extern template
class CaribouTopology<caribou::geometry::Quad<caribou::_2D, caribou::Linear>>;

// Quad 2D quadratic specialization
template<>
auto CaribouTopology<caribou::geometry::Quad<caribou::_2D, caribou::Quadratic>>::templateName(
        const CaribouTopology<caribou::geometry::Quad<caribou::_2D, caribou::Quadratic>> *) -> std::string;

extern template
class CaribouTopology<caribou::geometry::Quad<caribou::_2D, caribou::Quadratic>>;

// Quad 3D linear specialization
template<>
auto CaribouTopology<caribou::geometry::Quad<caribou::_3D, caribou::Linear>>::templateName(
        const CaribouTopology<caribou::geometry::Quad<caribou::_3D, caribou::Linear>> *) -> std::string;

extern template
class CaribouTopology<caribou::geometry::Quad<caribou::_3D, caribou::Linear>>;

// Quad 3D quadratic specialization
template<>
auto CaribouTopology<caribou::geometry::Quad<caribou::_3D, caribou::Quadratic>>::templateName(
        const CaribouTopology<caribou::geometry::Quad<caribou::_3D, caribou::Quadratic>> *) -> std::string;

extern template
class CaribouTopology<caribou::geometry::Quad<caribou::_3D, caribou::Quadratic>>;

} // namespace SofaCaribou::topology