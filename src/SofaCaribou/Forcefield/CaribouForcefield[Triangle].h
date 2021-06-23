#pragma once
#include <SofaCaribou/Forcefield/CaribouForcefield.h>
#include <SofaCaribou/Topology/CaribouTopology[Triangle].h>
#include <Caribou/Geometry/Triangle.h>

namespace SofaCaribou::forcefield {

// Triangle linear specialization
// 2D
template <> auto CaribouForcefield<caribou::geometry::Triangle < caribou::_2D, caribou::Linear>>::templateName(const CaribouForcefield<caribou::geometry::Triangle < caribou::_2D, caribou::Linear>> *) -> std::string;
template <> void CaribouForcefield<caribou::geometry::Triangle < caribou::_2D, caribou::Linear>>::triangulate_face(const caribou::geometry::Triangle < caribou::_2D, caribou::Linear> & e, const std::size_t & face_id, std::vector<sofa::type::Vector3> & triangles_nodes);
extern template class CaribouForcefield<caribou::geometry::Triangle < caribou::_2D, caribou::Linear>>;

// 3D
template <> auto CaribouForcefield<caribou::geometry::Triangle < caribou::_3D, caribou::Linear>>::templateName(const CaribouForcefield<caribou::geometry::Triangle < caribou::_3D, caribou::Linear>> *) -> std::string;
template <> void CaribouForcefield<caribou::geometry::Triangle < caribou::_3D, caribou::Linear>>::triangulate_face(const caribou::geometry::Triangle < caribou::_3D, caribou::Linear> & e, const std::size_t & face_id, std::vector<sofa::type::Vector3> & triangles_nodes);
extern template class CaribouForcefield<caribou::geometry::Triangle < caribou::_3D, caribou::Linear>>;

// Triangle quadratic specialization

// 2D
template <> auto CaribouForcefield<caribou::geometry::Triangle < caribou::_2D, caribou::Quadratic>>::templateName(const CaribouForcefield<caribou::geometry::Triangle < caribou::_2D, caribou::Quadratic>> *) -> std::string;
template <> void CaribouForcefield<caribou::geometry::Triangle < caribou::_2D, caribou::Quadratic>>::triangulate_face(const caribou::geometry::Triangle < caribou::_2D, caribou::Quadratic> & e, const std::size_t & face_id, std::vector<sofa::type::Vector3> & triangles_nodes);
extern template class CaribouForcefield<caribou::geometry::Triangle < caribou::_2D, caribou::Quadratic>>;

// 3D
template <> auto CaribouForcefield<caribou::geometry::Triangle < caribou::_3D, caribou::Quadratic>>::templateName(const CaribouForcefield<caribou::geometry::Triangle < caribou::_3D, caribou::Quadratic>> *) -> std::string;
template <> void CaribouForcefield<caribou::geometry::Triangle < caribou::_3D, caribou::Quadratic>>::triangulate_face(const caribou::geometry::Triangle < caribou::_3D, caribou::Quadratic> & e, const std::size_t & face_id, std::vector<sofa::type::Vector3> & triangles_nodes);
extern template class CaribouForcefield<caribou::geometry::Triangle < caribou::_3D, caribou::Quadratic>>;

}