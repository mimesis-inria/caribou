#pragma once
#include <SofaCaribou/Forcefield/CaribouForcefield.h>
#include <SofaCaribou/Topology/CaribouTopology[Tetrahedron10].h>
#include <Caribou/Geometry/Tetrahedron10.h>

namespace SofaCaribou::forcefield {

// Tetrahedron quadratic specialization
template <> auto CaribouForcefield<caribou::geometry::Tetrahedron10>::templateName(const CaribouForcefield<caribou::geometry::Tetrahedron10> *) -> std::string;
template <> void CaribouForcefield<caribou::geometry::Tetrahedron10>::triangulate_face(const caribou::geometry::Tetrahedron10 & e, const std::size_t & face_id, std::vector<sofa::type::Vec3> & triangles_nodes);
extern template class CaribouForcefield<caribou::geometry::Tetrahedron10>;

}