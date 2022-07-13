#pragma once

#include <SofaCaribou/FEniCS/Material/FEniCS_Material.h>
#include <SofaCaribou/Topology/CaribouTopology[Tetrahedron].h>

DISABLE_ALL_WARNINGS_BEGIN
#include <sofa/version.h>
#include <sofa/defaulttype/VecTypes.h>
DISABLE_ALL_WARNINGS_END

#include <Caribou/Geometry/Tetrahedron.h>

#include <ufcx.h>
namespace SofaCaribou::material {

// Hexahedron linear specialization
template <> auto FEniCS_Material<caribou::geometry::Tetrahedron, sofa::defaulttype::Vec3Types>::FEniCS_F() -> ufcx_integral*;
template <> auto FEniCS_Material<caribou::geometry::Tetrahedron, sofa::defaulttype::Vec3Types>::FEniCS_F_bc() -> ufcx_integral*;
template <> auto FEniCS_Material<caribou::geometry::Tetrahedron, sofa::defaulttype::Vec3Types>::FEniCS_J() -> ufcx_integral*;
template <> auto FEniCS_Material<caribou::geometry::Tetrahedron, sofa::defaulttype::Vec3Types>::FEniCS_Pi() -> ufcx_integral*;
extern template class FEniCS_Material<caribou::geometry::Tetrahedron, sofa::defaulttype::Vec3Types>;

} // namespace SofaCaribou::mapping
