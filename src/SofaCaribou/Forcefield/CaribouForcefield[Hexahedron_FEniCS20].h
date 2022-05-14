#pragma once
#include <SofaCaribou/Forcefield/CaribouForcefield.h>
#include <SofaCaribou/Topology/CaribouTopology[Hexahedron_FEniCS20].h>
#include <Caribou/Geometry/Hexahedron_FEniCS20.h>

namespace SofaCaribou::forcefield {

// Hexahedron_FEniCS quadratic specialization
template <> auto CaribouForcefield<caribou::geometry::Hexahedron_FEniCS20>::templateName(const CaribouForcefield<caribou::geometry::Hexahedron_FEniCS20> *) -> std::string;
template <> void CaribouForcefield<caribou::geometry::Hexahedron_FEniCS20>::triangulate_face(const caribou::geometry::Hexahedron_FEniCS20 & e, const std::size_t & face_id, std::vector<sofa::type::Vector3> & triangles_nodes);
extern template class CaribouForcefield<caribou::geometry::Hexahedron_FEniCS20>;

}
