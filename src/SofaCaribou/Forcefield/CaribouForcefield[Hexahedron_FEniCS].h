#pragma once
#include <SofaCaribou/Forcefield/CaribouForcefield.h>
#include <SofaCaribou/Topology/CaribouTopology[Hexahedron_FEniCS].h>
#include <Caribou/Geometry/Hexahedron_FEniCS.h>

namespace SofaCaribou::forcefield {

// Hexahedron_FEniCS linear specialization
template <> auto CaribouForcefield<caribou::geometry::Hexahedron_FEniCS < caribou::Linear>>::templateName(const CaribouForcefield<caribou::geometry::Hexahedron_FEniCS < caribou::Linear>> *) -> std::string;
template <> void CaribouForcefield<caribou::geometry::Hexahedron_FEniCS < caribou::Linear>>::triangulate_face(const caribou::geometry::Hexahedron_FEniCS < caribou::Linear> & e, const std::size_t & face_id, std::vector<sofa::type::Vector3> & triangles_nodes);
extern template class CaribouForcefield<caribou::geometry::Hexahedron_FEniCS < caribou::Linear>>;

// Hexahedron_FEniCS quadratic specialization
template <> auto CaribouForcefield<caribou::geometry::Hexahedron_FEniCS < caribou::Quadratic>>::templateName(const CaribouForcefield<caribou::geometry::Hexahedron_FEniCS < caribou::Quadratic>> *) -> std::string;
template <> void CaribouForcefield<caribou::geometry::Hexahedron_FEniCS < caribou::Quadratic>>::triangulate_face(const caribou::geometry::Hexahedron_FEniCS < caribou::Quadratic> & e, const std::size_t & face_id, std::vector<sofa::type::Vector3> & triangles_nodes);
extern template class CaribouForcefield<caribou::geometry::Hexahedron_FEniCS < caribou::Quadratic>>;

}
