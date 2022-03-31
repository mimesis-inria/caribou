#pragma once

#include <SofaCaribou/Topology/CaribouTopology.h>
#include <Caribou/Geometry/Hexahedron_FEniCS.h>

namespace SofaCaribou::topology {

// Hexahedron_FEniCS linear specialization
template<>
auto CaribouTopology<caribou::geometry::Hexahedron_FEniCS>::templateName(
        const CaribouTopology<caribou::geometry::Hexahedron_FEniCS> *) -> std::string;

template <> auto CaribouTopology<caribou::geometry::Hexahedron_FEniCS>::mesh_is_compatible(
        const sofa::core::topology::BaseMeshTopology * topology) -> bool;

template <> auto CaribouTopology<caribou::geometry::Hexahedron_FEniCS>::get_indices_data_from(
        const sofa::core::topology::BaseMeshTopology * topology) -> sofa::core::objectmodel::BaseData *;

extern template
class CaribouTopology<caribou::geometry::Hexahedron_FEniCS>;

} // namespace SofaCaribou::topology
