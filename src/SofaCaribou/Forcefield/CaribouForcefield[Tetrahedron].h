#pragma once
#include <SofaCaribou/Forcefield/CaribouForcefield.h>
#include <SofaCaribou/Topology/CaribouTopology[Tetrahedron].h>
#include <Caribou/Geometry/Tetrahedron.h>

namespace SofaCaribou::forcefield {

// Tetrahedron linear specialization
template <> auto CaribouForcefield<caribou::geometry::Tetrahedron>::templateName(const CaribouForcefield<caribou::geometry::Tetrahedron> *) -> std::string;
template <> void CaribouForcefield<caribou::geometry::Tetrahedron>::triangulate_face(const caribou::geometry::Tetrahedron & e, const std::size_t & face_id, std::vector<sofa::type::Vec3> & triangles_nodes);
extern template class CaribouForcefield<caribou::geometry::Tetrahedron>;

}