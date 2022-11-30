#pragma once

#include <SofaCaribou/Topology/CaribouTopology.h>
#include <Caribou/Geometry/Hexahedron.h>

namespace SofaCaribou::topology {

// Hexahedron linear specialization
template<>
auto CaribouTopology<caribou::geometry::Hexahedron>::GetCustomTemplateName() -> std::string;

template <> auto CaribouTopology<caribou::geometry::Hexahedron>::mesh_is_compatible(
        const sofa::core::topology::BaseMeshTopology * topology) -> bool;

template <> auto CaribouTopology<caribou::geometry::Hexahedron>::get_indices_data_from(
        const sofa::core::topology::BaseMeshTopology * topology) -> sofa::core::objectmodel::BaseData *;

extern template
class CaribouTopology<caribou::geometry::Hexahedron>;

} // namespace SofaCaribou::topology