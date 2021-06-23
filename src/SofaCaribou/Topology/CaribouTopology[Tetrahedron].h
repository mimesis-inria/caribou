#pragma once

#include <SofaCaribou/Topology/CaribouTopology.h>
#include <Caribou/Geometry/Tetrahedron.h>

namespace SofaCaribou::topology {

// Tetrahedron linear specialization
template<>
auto CaribouTopology<caribou::geometry::Tetrahedron<caribou::Linear>>::templateName(
        const CaribouTopology<caribou::geometry::Tetrahedron<caribou::Linear>> *) -> std::string;

template <> auto CaribouTopology<caribou::geometry::Tetrahedron < caribou::Linear>>::mesh_is_compatible(
        const sofa::core::topology::BaseMeshTopology * topology) -> bool;

template <> auto CaribouTopology<caribou::geometry::Tetrahedron < caribou::Linear>>::get_indices_data_from(
        const sofa::core::topology::BaseMeshTopology * topology) -> sofa::core::objectmodel::BaseData *;

extern template
class CaribouTopology<caribou::geometry::Tetrahedron<caribou::Linear>>;

// Tetrahedron quadratic specialization
template<>
auto CaribouTopology<caribou::geometry::Tetrahedron<caribou::Quadratic>>::templateName(
        const CaribouTopology<caribou::geometry::Tetrahedron<caribou::Quadratic>> *) -> std::string;

template <> auto CaribouTopology<caribou::geometry::Tetrahedron < caribou::Quadratic>>::mesh_is_compatible(
        const sofa::core::topology::BaseMeshTopology * topology) -> bool;

extern template
class CaribouTopology<caribou::geometry::Tetrahedron<caribou::Quadratic>>;

} // namespace SofaCaribou::topology