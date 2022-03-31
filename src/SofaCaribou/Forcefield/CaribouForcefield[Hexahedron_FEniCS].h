#pragma once
#include <SofaCaribou/Forcefield/CaribouForcefield.h>
#include <SofaCaribou/Topology/CaribouTopology[Hexahedron_FEniCS].h>
#include <Caribou/Geometry/Hexahedron_FEniCS.h>

namespace SofaCaribou::forcefield {

// Hexahedron_FEniCS linear specialization
template <> auto CaribouForcefield<caribou::geometry::Hexahedron_FEniCS>::templateName(const CaribouForcefield<caribou::geometry::Hexahedron_FEniCS> *) -> std::string;
template <> void CaribouForcefield<caribou::geometry::Hexahedron_FEniCS>::triangulate_face(const caribou::geometry::Hexahedron_FEniCS & e, const std::size_t & face_id, std::vector<sofa::type::Vector3> & triangles_nodes);
extern template class CaribouForcefield<caribou::geometry::Hexahedron_FEniCS>;

}
