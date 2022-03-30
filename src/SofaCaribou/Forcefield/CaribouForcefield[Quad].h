#pragma once
#include <SofaCaribou/Forcefield/CaribouForcefield.h>
#include <SofaCaribou/Topology/CaribouTopology[Quad].h>
#include <Caribou/Geometry/Quad.h>

namespace SofaCaribou::forcefield {

// Quad linear specialization
// 2D
template <> auto CaribouForcefield<caribou::geometry::Quad<caribou::_2D>>::templateName(const CaribouForcefield<caribou::geometry::Quad<caribou::_2D>> *) -> std::string;
template <> void CaribouForcefield<caribou::geometry::Quad<caribou::_2D>>::triangulate_face(const caribou::geometry::Quad<caribou::_2D> & e, const std::size_t & face_id, std::vector<sofa::type::Vector3> & triangles_nodes);
extern template class CaribouForcefield<caribou::geometry::Quad<caribou::_2D>>;

// 3D
template <> auto CaribouForcefield<caribou::geometry::Quad<caribou::_3D>>::templateName(const CaribouForcefield<caribou::geometry::Quad<caribou::_3D>> *) -> std::string;
template <> void CaribouForcefield<caribou::geometry::Quad<caribou::_3D>>::triangulate_face(const caribou::geometry::Quad<caribou::_3D> & e, const std::size_t & face_id, std::vector<sofa::type::Vector3> & triangles_nodes);
extern template class CaribouForcefield<caribou::geometry::Quad<caribou::_3D>>;

}