#pragma once

#include <Caribou/Geometry/Hexahedron.h>

#include <SofaCaribou/FEniCS/Material/NeoHookean/NeoHookeanMaterial_FEniCS.h>
#include <SofaCaribou/Topology/CaribouTopology[Tetrahedron10].h>

DISABLE_ALL_WARNINGS_BEGIN
#include <sofa/version.h>
#include <sofa/defaulttype/VecTypes.h>
DISABLE_ALL_WARNINGS_END

#include <Caribou/Geometry/Tetrahedron10.h>

#include <ufcx.h>
namespace SofaCaribou::material {

// Hexahedron linear specialization
template <> auto NeoHookeanMaterial_FEniCS<caribou::geometry::Tetrahedron10 , sofa::defaulttype::Vec3Types>::FEniCS_F() -> ufcx_integral*;
template <> auto NeoHookeanMaterial_FEniCS<caribou::geometry::Tetrahedron10 , sofa::defaulttype::Vec3Types>::FEniCS_J() -> ufcx_integral*;
template <> auto NeoHookeanMaterial_FEniCS<caribou::geometry::Tetrahedron10 , sofa::defaulttype::Vec3Types>::FEniCS_Pi() -> ufcx_integral*;
extern template class NeoHookeanMaterial_FEniCS<caribou::geometry::Tetrahedron10, sofa::defaulttype::Vec3Types>;

} // namespace SofaCaribou::mapping
