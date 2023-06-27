#pragma once

#include <SofaCaribou/Topology/CaribouTopology.h>
#include <Caribou/Geometry/Tetrahedron10.h>

namespace SofaCaribou::topology {

// Tetrahedron quadratic specialization
template<>
auto CaribouTopology<caribou::geometry::Tetrahedron10>::GetCustomTemplateName() -> std::string;

template <> auto CaribouTopology<caribou::geometry::Tetrahedron10>::mesh_is_compatible(
        const sofa::core::topology::BaseMeshTopology * topology) -> bool;

extern template
class CaribouTopology<caribou::geometry::Tetrahedron10>;

} // namespace SofaCaribou::topology