#pragma once

#include <SofaCaribou/Topology/CaribouTopology.h>
#include <Caribou/Geometry/Tetrahedron.h>

namespace SofaCaribou::topology {

// Tetrahedron linear specialization
template<>
auto CaribouTopology<caribou::geometry::Tetrahedron>::GetCustomTemplateName() -> std::string;

template <> auto CaribouTopology<caribou::geometry::Tetrahedron>::mesh_is_compatible(
        const sofa::core::topology::BaseMeshTopology * topology) -> bool;

template <> auto CaribouTopology<caribou::geometry::Tetrahedron>::get_indices_data_from(
        const sofa::core::topology::BaseMeshTopology * topology) -> sofa::core::objectmodel::BaseData *;

extern template
class CaribouTopology<caribou::geometry::Tetrahedron>;

} // namespace SofaCaribou::topology