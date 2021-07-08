#pragma once
#include <SofaCaribou/Forcefield/CaribouForcefield.h>
#include <SofaCaribou/Topology/CaribouTopology[Tetrahedron].h>
#include <Caribou/Geometry/Tetrahedron.h>

namespace SofaCaribou::forcefield {

// Tetrahedron linear specialization
template <> auto CaribouForcefield<caribou::geometry::Tetrahedron < caribou::Linear>>::get_indices_from(const sofa::core::topology::BaseMeshTopology * topology) -> sofa::core::objectmodel::BaseData *;
template <> auto CaribouForcefield<caribou::geometry::Tetrahedron < caribou::Linear>>::mesh_is_compatible(const sofa::core::topology::BaseMeshTopology * topology) -> bool;
template <> auto CaribouForcefield<caribou::geometry::Tetrahedron < caribou::Linear>>::templateName(const CaribouForcefield<caribou::geometry::Tetrahedron < caribou::Linear>> *) -> std::string;
template <> void CaribouForcefield<caribou::geometry::Tetrahedron < caribou::Linear>>::triangulate_face(const caribou::geometry::Tetrahedron < caribou::Linear> & e, const std::size_t & face_id, std::vector<sofa::type::Vector3> & triangles_nodes);
extern template class CaribouForcefield<caribou::geometry::Tetrahedron < caribou::Linear>>;

// Tetrahedron quadratic specialization
template <> auto CaribouForcefield<caribou::geometry::Tetrahedron < caribou::Quadratic>>::get_indices_from(const sofa::core::topology::BaseMeshTopology * topology) -> sofa::core::objectmodel::BaseData *;
template <> auto CaribouForcefield<caribou::geometry::Tetrahedron < caribou::Quadratic>>::mesh_is_compatible(const sofa::core::topology::BaseMeshTopology * topology) -> bool;
template <> auto CaribouForcefield<caribou::geometry::Tetrahedron < caribou::Quadratic>>::templateName(const CaribouForcefield<caribou::geometry::Tetrahedron < caribou::Quadratic>> *) -> std::string;
template <> void CaribouForcefield<caribou::geometry::Tetrahedron < caribou::Quadratic>>::triangulate_face(const caribou::geometry::Tetrahedron < caribou::Quadratic> & e, const std::size_t & face_id, std::vector<sofa::type::Vector3> & triangles_nodes);
extern template class CaribouForcefield<caribou::geometry::Tetrahedron < caribou::Quadratic>>;

}