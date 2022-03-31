#pragma once

#include <Caribou/Geometry/Hexahedron_FEniCS.h>

#include <SofaCaribou/FEniCS/Material/SaintVenantKirchhoff/SaintVenantKirchhoffMaterial_FEniCS.h>
#include <SofaCaribou/Topology/CaribouTopology[Hexahedron_FEniCS].h>

DISABLE_ALL_WARNINGS_BEGIN
#include <sofa/version.h>
#include <sofa/defaulttype/VecTypes.h>
DISABLE_ALL_WARNINGS_END

#include <ufcx.h>
#include <Caribou/Geometry/Hexahedron_FEniCS.h>

namespace SofaCaribou::material {

// Hexahedron_FEniCS linear specialization
template <> auto SaintVenantKirchhoffMaterial_FEniCS<caribou::geometry::Hexahedron_FEniCS, sofa::defaulttype::Vec3Types>::FEniCS_F() -> ufcx_integral*;
template <> auto SaintVenantKirchhoffMaterial_FEniCS<caribou::geometry::Hexahedron_FEniCS, sofa::defaulttype::Vec3Types>::FEniCS_J() -> ufcx_integral*;
extern template class SaintVenantKirchhoffMaterial_FEniCS<caribou::geometry::Hexahedron_FEniCS, sofa::defaulttype::Vec3Types>;

} // namespace SofaCaribou::mapping
