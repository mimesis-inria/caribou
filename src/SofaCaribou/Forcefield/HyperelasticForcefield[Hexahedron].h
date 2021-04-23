#pragma once
#include <SofaCaribou/Forcefield/HyperelasticForcefield.h>
#include <Caribou/Geometry/Hexahedron.h>

namespace SofaCaribou::forcefield {

// Hexahedron linear specialization
template <> void HyperelasticForcefield<caribou::geometry::Hexahedron < caribou::Linear>>::create_domain_from(sofa::core::topology::BaseMeshTopology * topology) ;
template <> auto HyperelasticForcefield<caribou::geometry::Hexahedron < caribou::Linear>>::mesh_is_compatible(const sofa::core::topology::BaseMeshTopology * topology) -> bool;
template <> auto HyperelasticForcefield<caribou::geometry::Hexahedron < caribou::Linear>>::templateName(const HyperelasticForcefield<caribou::geometry::Hexahedron < caribou::Linear>> *) -> std::string;
extern template class HyperelasticForcefield<caribou::geometry::Hexahedron < caribou::Linear>>;

// Hexahedron quadratic specialization
template <> void HyperelasticForcefield<caribou::geometry::Hexahedron < caribou::Quadratic>>::create_domain_from(sofa::core::topology::BaseMeshTopology * topology) ;
template <> auto HyperelasticForcefield<caribou::geometry::Hexahedron < caribou::Quadratic>>::mesh_is_compatible(const sofa::core::topology::BaseMeshTopology * topology) -> bool;
template <> auto HyperelasticForcefield<caribou::geometry::Hexahedron < caribou::Quadratic>>::templateName(const HyperelasticForcefield<caribou::geometry::Hexahedron < caribou::Quadratic>> *) -> std::string;
template <> void HyperelasticForcefield<caribou::geometry::Hexahedron < caribou::Quadratic>>::draw(const sofa::core::visual::VisualParams * vparams);
extern template class HyperelasticForcefield<caribou::geometry::Hexahedron < caribou::Quadratic>>;

}