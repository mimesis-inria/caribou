#pragma once
#include <SofaCaribou/Forcefield/CaribouForcefield.h>
#include <SofaCaribou/Topology/CaribouTopology[Hexahedron].h>
#include <Caribou/Geometry/Hexahedron.h>

namespace SofaCaribou::forcefield {

// Hexahedron linear specialization
template <> auto CaribouForcefield<caribou::geometry::Hexahedron>::templateName(const CaribouForcefield<caribou::geometry::Hexahedron> *) -> std::string;
template <> void CaribouForcefield<caribou::geometry::Hexahedron>::triangulate_face(const caribou::geometry::Hexahedron & e, const std::size_t & face_id, std::vector<sofa::type::Vec3> & triangles_nodes);
extern template class CaribouForcefield<caribou::geometry::Hexahedron>;

}