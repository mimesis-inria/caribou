#pragma once

#include <SofaCaribou/Topology/CaribouTopology.h>
#include <Caribou/Geometry/Hexahedron_FEniCS20.h>

namespace SofaCaribou::topology {

// Hexahedron_FEniCS quadratic specialization
template<>
auto CaribouTopology<caribou::geometry::Hexahedron_FEniCS20>::templateName(
        const CaribouTopology<caribou::geometry::Hexahedron_FEniCS20> *) -> std::string;

template <> auto CaribouTopology<caribou::geometry::Hexahedron_FEniCS20>::mesh_is_compatible(
        const sofa::core::topology::BaseMeshTopology * topology) -> bool;

template <> auto CaribouTopology<caribou::geometry::Hexahedron_FEniCS20>::get_indices_data_from(
        const sofa::core::topology::BaseMeshTopology * topology) -> sofa::core::objectmodel::BaseData *;

extern template
class CaribouTopology<caribou::geometry::Hexahedron_FEniCS20>;

} // namespace SofaCaribou::topology
