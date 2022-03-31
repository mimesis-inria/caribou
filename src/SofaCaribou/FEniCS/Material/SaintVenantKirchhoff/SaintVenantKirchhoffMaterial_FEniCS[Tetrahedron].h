#pragma once

#include <Caribou/Geometry/Hexahedron.h>

#include <SofaCaribou/FEniCS/Material/SaintVenantKirchhoff/SaintVenantKirchhoffMaterial_FEniCS.h>
#include <SofaCaribou/Topology/CaribouTopology[Tetrahedron].h>

DISABLE_ALL_WARNINGS_BEGIN
#include <sofa/version.h>
#include <sofa/defaulttype/VecTypes.h>
DISABLE_ALL_WARNINGS_END

#include <Caribou/Geometry/Tetrahedron.h>

#include <ufcx.h>
namespace SofaCaribou::material {

// Hexahedron linear specialization
template <> auto SaintVenantKirchhoffMaterial_FEniCS<caribou::geometry::Tetrahedron, sofa::defaulttype::Vec3Types>::FEniCS_F() -> ufcx_integral*;
template <> auto SaintVenantKirchhoffMaterial_FEniCS<caribou::geometry::Tetrahedron, sofa::defaulttype::Vec3Types>::FEniCS_J() -> ufcx_integral*;
extern template class SaintVenantKirchhoffMaterial_FEniCS<caribou::geometry::Tetrahedron, sofa::defaulttype::Vec3Types>;


} // namespace SofaCaribou::mapping
