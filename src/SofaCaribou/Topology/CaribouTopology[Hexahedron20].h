#pragma once

#include <SofaCaribou/Topology/CaribouTopology.h>
#include <Caribou/Geometry/Hexahedron20.h>

namespace SofaCaribou::topology {

// Hexahedron quadratic specialization
template<>
auto CaribouTopology<caribou::geometry::Hexahedron20>::templateName(
        const CaribouTopology<caribou::geometry::Hexahedron20> *) -> std::string;

template <> auto CaribouTopology<caribou::geometry::Hexahedron20>::mesh_is_compatible(
        const sofa::core::topology::BaseMeshTopology * topology) -> bool;

extern template
class CaribouTopology<caribou::geometry::Hexahedron20>;

} // namespace SofaCaribou::topology