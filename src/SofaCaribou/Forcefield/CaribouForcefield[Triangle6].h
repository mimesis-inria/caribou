#pragma once
#include <SofaCaribou/Forcefield/CaribouForcefield.h>
#include <SofaCaribou/Topology/CaribouTopology[Triangle6].h>
#include <Caribou/Geometry/Triangle6.h>

namespace SofaCaribou::forcefield {

// Triangle quadratic specialization

// 2D
template <> auto CaribouForcefield<caribou::geometry::Triangle<caribou::_2D>>::templateName(const CaribouForcefield<caribou::geometry::Triangle<caribou::_2D>> *) -> std::string;
template <> void CaribouForcefield<caribou::geometry::Triangle<caribou::_2D>>::triangulate_face(const caribou::geometry::Triangle<caribou::_2D> & e, const std::size_t & face_id, std::vector<sofa::type::Vec3> & triangles_nodes);
extern template class CaribouForcefield<caribou::geometry::Triangle<caribou::_2D>>;

// 3D
template <> auto CaribouForcefield<caribou::geometry::Triangle<caribou::_3D>>::templateName(const CaribouForcefield<caribou::geometry::Triangle<caribou::_3D>> *) -> std::string;
template <> void CaribouForcefield<caribou::geometry::Triangle<caribou::_3D>>::triangulate_face(const caribou::geometry::Triangle<caribou::_3D> & e, const std::size_t & face_id, std::vector<sofa::type::Vec3> & triangles_nodes);
extern template class CaribouForcefield<caribou::geometry::Triangle<caribou::_3D>>;

}