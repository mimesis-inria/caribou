#pragma once

#include <SofaCaribou/Topology/CaribouTopology.h>
#include <Caribou/Geometry/Hexahedron.h>

namespace SofaCaribou::topology {

// Hexahedron linear specialization
template<>
auto CaribouTopology<caribou::geometry::Hexahedron<caribou::Linear>>::templateName(
        const CaribouTopology<caribou::geometry::Hexahedron<caribou::Linear>> *) -> std::string;

extern template
class CaribouTopology<caribou::geometry::Hexahedron<caribou::Linear>>;

// Hexahedron quadratic specialization
template<>
auto CaribouTopology<caribou::geometry::Hexahedron<caribou::Quadratic>>::templateName(
        const CaribouTopology<caribou::geometry::Hexahedron<caribou::Quadratic>> *) -> std::string;

extern template
class CaribouTopology<caribou::geometry::Hexahedron<caribou::Quadratic>>;

} // namespace SofaCaribou::topology