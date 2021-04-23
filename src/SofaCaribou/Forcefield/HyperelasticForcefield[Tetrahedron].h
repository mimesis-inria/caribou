#pragma once
#include <SofaCaribou/Forcefield/HyperelasticForcefield.h>
#include <Caribou/Geometry/Tetrahedron.h>

namespace SofaCaribou::forcefield {

// Tetrahedron linear specialization
template <> void HyperelasticForcefield<caribou::geometry::Tetrahedron < caribou::Linear>>::create_domain_from(sofa::core::topology::BaseMeshTopology * topology) ;
template <> auto HyperelasticForcefield<caribou::geometry::Tetrahedron < caribou::Linear>>::mesh_is_compatible(const sofa::core::topology::BaseMeshTopology * topology) -> bool;
template <> auto HyperelasticForcefield<caribou::geometry::Tetrahedron < caribou::Linear>>::templateName(const HyperelasticForcefield<caribou::geometry::Tetrahedron < caribou::Linear>> *) -> std::string;
extern template class HyperelasticForcefield<caribou::geometry::Tetrahedron < caribou::Linear>>;

// Tetrahedron quadratic specialization
template <> void HyperelasticForcefield<caribou::geometry::Tetrahedron < caribou::Quadratic>>::create_domain_from(sofa::core::topology::BaseMeshTopology * topology) ;
template <> auto HyperelasticForcefield<caribou::geometry::Tetrahedron < caribou::Quadratic>>::mesh_is_compatible(const sofa::core::topology::BaseMeshTopology * topology) -> bool;
template <> auto HyperelasticForcefield<caribou::geometry::Tetrahedron < caribou::Quadratic>>::templateName(const HyperelasticForcefield<caribou::geometry::Tetrahedron < caribou::Quadratic>> *) -> std::string;
template <> void HyperelasticForcefield<caribou::geometry::Tetrahedron < caribou::Quadratic>>::draw(const sofa::core::visual::VisualParams * vparams);
extern template class HyperelasticForcefield<caribou::geometry::Tetrahedron < caribou::Quadratic>>;

}