#pragma once

#include <SofaCaribou/Topology/CaribouTopology.h>
#include <Caribou/Geometry/Triangle.h>

namespace SofaCaribou::topology {

// Triangle 2D linear specialization
template<>
auto CaribouTopology<caribou::geometry::Triangle<caribou::_2D>>::templateName(
        const CaribouTopology<caribou::geometry::Triangle<caribou::_2D>> *) -> std::string;

template <> auto CaribouTopology<caribou::geometry::Triangle<caribou::_2D>>::mesh_is_compatible(
        const sofa::core::topology::BaseMeshTopology * topology) -> bool;

template <> auto CaribouTopology<caribou::geometry::Triangle<caribou::_2D>>::get_indices_data_from(
        const sofa::core::topology::BaseMeshTopology * topology) -> sofa::core::objectmodel::BaseData *;

extern template
class CaribouTopology<caribou::geometry::Triangle<caribou::_2D>>;

// Triangle 3D linear specialization
template<>
auto CaribouTopology<caribou::geometry::Triangle<caribou::_3D>>::templateName(
        const CaribouTopology<caribou::geometry::Triangle<caribou::_3D>> *) -> std::string;

template <> auto CaribouTopology<caribou::geometry::Triangle<caribou::_3D>>::mesh_is_compatible(
        const sofa::core::topology::BaseMeshTopology * topology) -> bool;

template <> auto CaribouTopology<caribou::geometry::Triangle<caribou::_3D>>::get_indices_data_from(
        const sofa::core::topology::BaseMeshTopology * topology) -> sofa::core::objectmodel::BaseData *;

extern template
class CaribouTopology<caribou::geometry::Triangle<caribou::_3D>>;

} // namespace SofaCaribou::topology