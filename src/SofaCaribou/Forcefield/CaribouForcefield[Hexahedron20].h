#pragma once
#include <SofaCaribou/Forcefield/CaribouForcefield.h>
#include <SofaCaribou/Topology/CaribouTopology[Hexahedron20].h>
#include <Caribou/Geometry/Hexahedron20.h>

namespace SofaCaribou::forcefield {

// Hexahedron quadratic specialization
template <> auto CaribouForcefield<caribou::geometry::Hexahedron20>::templateName(const CaribouForcefield<caribou::geometry::Hexahedron20> *) -> std::string;
template <> void CaribouForcefield<caribou::geometry::Hexahedron20>::triangulate_face(const caribou::geometry::Hexahedron20 & e, const std::size_t & face_id, std::vector<sofa::type::Vec3> & triangles_nodes);
extern template class CaribouForcefield<caribou::geometry::Hexahedron20>;

}