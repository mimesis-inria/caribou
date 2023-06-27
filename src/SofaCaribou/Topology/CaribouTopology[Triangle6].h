#pragma once

#include <SofaCaribou/Topology/CaribouTopology.h>
#include <Caribou/Geometry/Triangle6.h>

namespace SofaCaribou::topology {

// Triangle 2D quadratic specialization
template<> auto CaribouTopology<caribou::geometry::Triangle6<caribou::_2D>>::GetCustomTemplateName() -> std::string;

template <> auto CaribouTopology<caribou::geometry::Triangle6<caribou::_2D>>::mesh_is_compatible(
        const sofa::core::topology::BaseMeshTopology * topology) -> bool;

extern template
class CaribouTopology<caribou::geometry::Triangle6<caribou::_2D>>;

// Triangle 3D quadratic specialization
template<> auto CaribouTopology<caribou::geometry::Triangle6<caribou::_3D>>::GetCustomTemplateName() -> std::string;

template <> auto CaribouTopology<caribou::geometry::Triangle6<caribou::_3D>>::mesh_is_compatible(
        const sofa::core::topology::BaseMeshTopology * topology) -> bool;

extern template
class CaribouTopology<caribou::geometry::Triangle6<caribou::_3D>>;

} // namespace SofaCaribou::topology