#pragma once

#include <SofaCaribou/Topology/CaribouTopology.h>
#include <Caribou/Geometry/Tetrahedron.h>

namespace SofaCaribou::topology {

// Tetrahedron linear specialization
template<>
auto CaribouTopology<caribou::geometry::Tetrahedron<caribou::Linear>>::templateName(
        const CaribouTopology<caribou::geometry::Tetrahedron<caribou::Linear>> *) -> std::string;

extern template
class CaribouTopology<caribou::geometry::Tetrahedron<caribou::Linear>>;

// Tetrahedron quadratic specialization
template<>
auto CaribouTopology<caribou::geometry::Tetrahedron<caribou::Quadratic>>::templateName(
        const CaribouTopology<caribou::geometry::Tetrahedron<caribou::Quadratic>> *) -> std::string;

extern template
class CaribouTopology<caribou::geometry::Tetrahedron<caribou::Quadratic>>;

} // namespace SofaCaribou::topology