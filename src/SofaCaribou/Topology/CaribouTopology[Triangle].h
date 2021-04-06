#pragma once

#include <SofaCaribou/Topology/CaribouTopology.h>
#include <Caribou/Geometry/Triangle.h>

namespace SofaCaribou::topology {

// Triangle 2D linear specialization
template<>
auto CaribouTopology<caribou::geometry::Triangle<caribou::_2D, caribou::Linear>>::templateName(
        const CaribouTopology<caribou::geometry::Triangle<caribou::_2D, caribou::Linear>> *) -> std::string;

extern template
class CaribouTopology<caribou::geometry::Triangle<caribou::_2D, caribou::Linear>>;

// Triangle 2D quadratic specialization
template<>
auto CaribouTopology<caribou::geometry::Triangle<caribou::_2D, caribou::Quadratic>>::templateName(
        const CaribouTopology<caribou::geometry::Triangle<caribou::_2D, caribou::Quadratic>> *) -> std::string;

extern template
class CaribouTopology<caribou::geometry::Triangle<caribou::_2D, caribou::Quadratic>>;

// Triangle 3D linear specialization
template<>
auto CaribouTopology<caribou::geometry::Triangle<caribou::_3D, caribou::Linear>>::templateName(
        const CaribouTopology<caribou::geometry::Triangle<caribou::_3D, caribou::Linear>> *) -> std::string;

extern template
class CaribouTopology<caribou::geometry::Triangle<caribou::_3D, caribou::Linear>>;

// Triangle 3D quadratic specialization
template<>
auto CaribouTopology<caribou::geometry::Triangle<caribou::_3D, caribou::Quadratic>>::templateName(
        const CaribouTopology<caribou::geometry::Triangle<caribou::_3D, caribou::Quadratic>> *) -> std::string;

extern template
class CaribouTopology<caribou::geometry::Triangle<caribou::_3D, caribou::Quadratic>>;

} // namespace SofaCaribou::topology