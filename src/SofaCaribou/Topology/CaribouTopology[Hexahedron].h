#pragma once

#include <SofaCaribou/Topology/CaribouTopology.h>
#include <Caribou/Geometry/Hexahedron.h>

namespace SofaCaribou::topology {

// Hexahedron linear specialization
template<>
auto CaribouTopology<caribou::geometry::Hexahedron<caribou::Linear>>::templateName(
        const CaribouTopology<caribou::geometry::Hexahedron<caribou::Linear>> *) -> std::string;

template <> auto CaribouTopology<caribou::geometry::Hexahedron < caribou::Linear>>::mesh_is_compatible(
        const sofa::core::topology::BaseMeshTopology * topology) -> bool;

template <> auto CaribouTopology<caribou::geometry::Hexahedron < caribou::Linear>>::get_indices_data_from(
        const sofa::core::topology::BaseMeshTopology * topology) -> sofa::core::objectmodel::BaseData *;

extern template
class CaribouTopology<caribou::geometry::Hexahedron<caribou::Linear>>;

// Hexahedron quadratic specialization
template<>
auto CaribouTopology<caribou::geometry::Hexahedron<caribou::Quadratic>>::templateName(
        const CaribouTopology<caribou::geometry::Hexahedron<caribou::Quadratic>> *) -> std::string;

template <> auto CaribouTopology<caribou::geometry::Hexahedron < caribou::Quadratic>>::mesh_is_compatible(
        const sofa::core::topology::BaseMeshTopology * topology) -> bool;

extern template
class CaribouTopology<caribou::geometry::Hexahedron<caribou::Quadratic>>;

} // namespace SofaCaribou::topology