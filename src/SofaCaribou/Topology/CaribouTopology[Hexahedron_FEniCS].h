#pragma once

#include <SofaCaribou/Topology/CaribouTopology.h>
#include <Caribou/Geometry/Hexahedron_FEniCS.h>

namespace SofaCaribou::topology {

// Hexahedron_FEniCS linear specialization
template<>
auto CaribouTopology<caribou::geometry::Hexahedron_FEniCS<caribou::Linear>>::templateName(
        const CaribouTopology<caribou::geometry::Hexahedron_FEniCS<caribou::Linear>> *) -> std::string;

template <> auto CaribouTopology<caribou::geometry::Hexahedron_FEniCS < caribou::Linear>>::mesh_is_compatible(
        const sofa::core::topology::BaseMeshTopology * topology) -> bool;

template <> auto CaribouTopology<caribou::geometry::Hexahedron_FEniCS < caribou::Linear>>::get_indices_data_from(
        const sofa::core::topology::BaseMeshTopology * topology) -> sofa::core::objectmodel::BaseData *;

extern template
class CaribouTopology<caribou::geometry::Hexahedron_FEniCS<caribou::Linear>>;

// Hexahedron_FEniCS quadratic specialization
template<>
auto CaribouTopology<caribou::geometry::Hexahedron_FEniCS<caribou::Quadratic>>::templateName(
        const CaribouTopology<caribou::geometry::Hexahedron_FEniCS<caribou::Quadratic>> *) -> std::string;

template <> auto CaribouTopology<caribou::geometry::Hexahedron_FEniCS < caribou::Quadratic>>::mesh_is_compatible(
        const sofa::core::topology::BaseMeshTopology * topology) -> bool;

extern template
class CaribouTopology<caribou::geometry::Hexahedron_FEniCS<caribou::Quadratic>>;

} // namespace SofaCaribou::topology
