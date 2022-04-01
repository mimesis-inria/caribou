#pragma once

#include <Caribou/Geometry/Hexahedron.h>

#include <SofaCaribou/FEniCS/Material/FEniCS_Material.h>
#include <SofaCaribou/Topology/CaribouTopology[Hexahedron].h>

DISABLE_ALL_WARNINGS_BEGIN
#include <sofa/version.h>
#include <sofa/defaulttype/VecTypes.h>
DISABLE_ALL_WARNINGS_END

#include <ufcx.h>
#include <Caribou/Geometry/Hexahedron.h>

namespace SofaCaribou::material {

// Hexahedron linear specialization
template <> auto FEniCS_Material<caribou::geometry::Hexahedron <caribou::Linear>, sofa::defaulttype::Vec3Types>::FEniCS_F() -> ufcx_integral*;
template <> auto FEniCS_Material<caribou::geometry::Hexahedron <caribou::Linear>, sofa::defaulttype::Vec3Types>::FEniCS_J() -> ufcx_integral*;
extern template class FEniCS_Material<caribou::geometry::Hexahedron<caribou::Linear>, sofa::defaulttype::Vec3Types>;

template <> auto FEniCS_Material<caribou::geometry::Hexahedron <caribou::Quadratic>, sofa::defaulttype::Vec3Types>::FEniCS_F() -> ufcx_integral*;
template <> auto FEniCS_Material<caribou::geometry::Hexahedron <caribou::Quadratic>, sofa::defaulttype::Vec3Types>::FEniCS_J() -> ufcx_integral*;
extern template class FEniCS_Material<caribou::geometry::Hexahedron<caribou::Quadratic>, sofa::defaulttype::Vec3Types>;


} // namespace SofaCaribou::mapping
