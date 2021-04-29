#pragma once
#include <SofaCaribou/Forcefield/CaribouForcefield.h>
#include <SofaCaribou/Topology/CaribouTopology[Quad].h>
#include <Caribou/Geometry/Quad.h>

namespace SofaCaribou::forcefield {

// Quad linear specialization
// 2D
template <> auto CaribouForcefield<caribou::geometry::Quad < caribou::_2D, caribou::Linear>>::get_indices_from(const sofa::core::topology::BaseMeshTopology * topology) -> sofa::core::objectmodel::BaseData *;
template <> auto CaribouForcefield<caribou::geometry::Quad < caribou::_2D, caribou::Linear>>::mesh_is_compatible(const sofa::core::topology::BaseMeshTopology * topology) -> bool;
template <> auto CaribouForcefield<caribou::geometry::Quad < caribou::_2D, caribou::Linear>>::templateName(const CaribouForcefield<caribou::geometry::Quad < caribou::_2D, caribou::Linear>> *) -> std::string;
template <> void CaribouForcefield<caribou::geometry::Quad < caribou::_2D, caribou::Linear>>::triangulate_face(const caribou::geometry::Quad < caribou::_2D, caribou::Linear> & e, const std::size_t & face_id, std::vector<sofa::defaulttype::Vector3> & triangles_nodes);
extern template class CaribouForcefield<caribou::geometry::Quad < caribou::_2D, caribou::Linear>>;

// 3D
template <> auto CaribouForcefield<caribou::geometry::Quad < caribou::_3D, caribou::Linear>>::get_indices_from(const sofa::core::topology::BaseMeshTopology * topology) -> sofa::core::objectmodel::BaseData *;
template <> auto CaribouForcefield<caribou::geometry::Quad < caribou::_3D, caribou::Linear>>::mesh_is_compatible(const sofa::core::topology::BaseMeshTopology * topology) -> bool;
template <> auto CaribouForcefield<caribou::geometry::Quad < caribou::_3D, caribou::Linear>>::templateName(const CaribouForcefield<caribou::geometry::Quad < caribou::_3D, caribou::Linear>> *) -> std::string;
template <> void CaribouForcefield<caribou::geometry::Quad < caribou::_3D, caribou::Linear>>::triangulate_face(const caribou::geometry::Quad < caribou::_3D, caribou::Linear> & e, const std::size_t & face_id, std::vector<sofa::defaulttype::Vector3> & triangles_nodes);
extern template class CaribouForcefield<caribou::geometry::Quad < caribou::_3D, caribou::Linear>>;

// Quad quadratic specialization

// 2D
template <> auto CaribouForcefield<caribou::geometry::Quad < caribou::_2D, caribou::Quadratic>>::get_indices_from(const sofa::core::topology::BaseMeshTopology * topology) -> sofa::core::objectmodel::BaseData *;
template <> auto CaribouForcefield<caribou::geometry::Quad < caribou::_2D, caribou::Quadratic>>::mesh_is_compatible(const sofa::core::topology::BaseMeshTopology * topology) -> bool;
template <> auto CaribouForcefield<caribou::geometry::Quad < caribou::_2D, caribou::Quadratic>>::templateName(const CaribouForcefield<caribou::geometry::Quad < caribou::_2D, caribou::Quadratic>> *) -> std::string;
template <> void CaribouForcefield<caribou::geometry::Quad < caribou::_2D, caribou::Quadratic>>::triangulate_face(const caribou::geometry::Quad < caribou::_2D, caribou::Quadratic> & e, const std::size_t & face_id, std::vector<sofa::defaulttype::Vector3> & triangles_nodes);
extern template class CaribouForcefield<caribou::geometry::Quad < caribou::_2D, caribou::Quadratic>>;

// 3D
template <> auto CaribouForcefield<caribou::geometry::Quad < caribou::_3D, caribou::Quadratic>>::get_indices_from(const sofa::core::topology::BaseMeshTopology * topology) -> sofa::core::objectmodel::BaseData *;
template <> auto CaribouForcefield<caribou::geometry::Quad < caribou::_3D, caribou::Quadratic>>::mesh_is_compatible(const sofa::core::topology::BaseMeshTopology * topology) -> bool;
template <> auto CaribouForcefield<caribou::geometry::Quad < caribou::_3D, caribou::Quadratic>>::templateName(const CaribouForcefield<caribou::geometry::Quad < caribou::_3D, caribou::Quadratic>> *) -> std::string;
template <> void CaribouForcefield<caribou::geometry::Quad < caribou::_3D, caribou::Quadratic>>::triangulate_face(const caribou::geometry::Quad < caribou::_3D, caribou::Quadratic> & e, const std::size_t & face_id, std::vector<sofa::defaulttype::Vector3> & triangles_nodes);
extern template class CaribouForcefield<caribou::geometry::Quad < caribou::_3D, caribou::Quadratic>>;

}