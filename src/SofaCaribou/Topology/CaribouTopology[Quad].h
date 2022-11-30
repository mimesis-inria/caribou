#pragma once

#include <SofaCaribou/Topology/CaribouTopology.h>
#include <Caribou/Geometry/Quad.h>

namespace SofaCaribou::topology {

// Quad 2D linear specialization
template<>
auto CaribouTopology<caribou::geometry::Quad<caribou::_2D>>::GetCustomTemplateName() -> std::string;

template <> auto CaribouTopology<caribou::geometry::Quad<caribou::_2D>>::mesh_is_compatible(
        const sofa::core::topology::BaseMeshTopology * topology) -> bool;

template <> auto CaribouTopology<caribou::geometry::Quad<caribou::_2D>>::get_indices_data_from(
        const sofa::core::topology::BaseMeshTopology * topology) -> sofa::core::objectmodel::BaseData *;

extern template
class CaribouTopology<caribou::geometry::Quad<caribou::_2D>>;

// Quad 3D linear specialization
template<>
auto CaribouTopology<caribou::geometry::Quad<caribou::_3D>>::GetCustomTemplateName() -> std::string;

template <> auto CaribouTopology<caribou::geometry::Quad <caribou::_3D>>::mesh_is_compatible(
        const sofa::core::topology::BaseMeshTopology * topology) -> bool;

template <> auto CaribouTopology<caribou::geometry::Quad <caribou::_3D>>::get_indices_data_from(
        const sofa::core::topology::BaseMeshTopology * topology) -> sofa::core::objectmodel::BaseData *;

extern template
class CaribouTopology<caribou::geometry::Quad<caribou::_3D>>;

} // namespace SofaCaribou::topology