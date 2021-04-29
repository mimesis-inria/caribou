#pragma once
#include <SofaCaribou/Forcefield/CaribouForcefield.h>
#include <SofaCaribou/Topology/CaribouTopology[Hexahedron].h>
#include <Caribou/Geometry/Hexahedron.h>

namespace SofaCaribou::forcefield {

// Hexahedron linear specialization
template <> auto CaribouForcefield<caribou::geometry::Hexahedron < caribou::Linear>>::get_indices_from(const sofa::core::topology::BaseMeshTopology * topology) -> sofa::core::objectmodel::BaseData *;
template <> auto CaribouForcefield<caribou::geometry::Hexahedron < caribou::Linear>>::mesh_is_compatible(const sofa::core::topology::BaseMeshTopology * topology) -> bool;
template <> auto CaribouForcefield<caribou::geometry::Hexahedron < caribou::Linear>>::templateName(const CaribouForcefield<caribou::geometry::Hexahedron < caribou::Linear>> *) -> std::string;
template <> void CaribouForcefield<caribou::geometry::Hexahedron < caribou::Linear>>::triangulate_face(const caribou::geometry::Hexahedron < caribou::Linear> & e, const std::size_t & face_id, std::vector<sofa::defaulttype::Vector3> & triangles_nodes);
extern template class CaribouForcefield<caribou::geometry::Hexahedron < caribou::Linear>>;

// Hexahedron quadratic specialization
template <> auto CaribouForcefield<caribou::geometry::Hexahedron < caribou::Quadratic>>::get_indices_from(const sofa::core::topology::BaseMeshTopology * topology) -> sofa::core::objectmodel::BaseData *;
template <> auto CaribouForcefield<caribou::geometry::Hexahedron < caribou::Quadratic>>::mesh_is_compatible(const sofa::core::topology::BaseMeshTopology * topology) -> bool;
template <> auto CaribouForcefield<caribou::geometry::Hexahedron < caribou::Quadratic>>::templateName(const CaribouForcefield<caribou::geometry::Hexahedron < caribou::Quadratic>> *) -> std::string;
template <> void CaribouForcefield<caribou::geometry::Hexahedron < caribou::Quadratic>>::triangulate_face(const caribou::geometry::Hexahedron < caribou::Quadratic> & e, const std::size_t & face_id, std::vector<sofa::defaulttype::Vector3> & triangles_nodes);
extern template class CaribouForcefield<caribou::geometry::Hexahedron < caribou::Quadratic>>;

}