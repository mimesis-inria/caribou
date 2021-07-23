#include "HyperelasticForcefieldRecomputeF.inl"
#include "HyperelasticForcefieldStoreF.inl"
#include "HyperelasticForcefieldStoreFAndS.inl"

#include <SofaCaribou/Forcefield/CaribouForcefield[Hexahedron].h>
#include <SofaCaribou/Forcefield/CaribouForcefield[Tetrahedron].h>
#include <Caribou/Geometry/Hexahedron.h>
#include <Caribou/Geometry/Tetrahedron.h>

DISABLE_ALL_WARNINGS_BEGIN
#include <sofa/core/ObjectFactory.h>
DISABLE_ALL_WARNINGS_END

using sofa::core::RegisterObject;
using namespace caribou::geometry;
using namespace caribou;
using namespace SofaCaribou::benchmark::forcefield;

[[maybe_unused]]
static int _c1_ = RegisterObject("Caribou hyperelastic force field")
.add<HyperelasticForcefieldRecomputeF<Hexahedron<Linear>>>()
.add<HyperelasticForcefieldRecomputeF<Tetrahedron<Linear>>>();

static int _c2_ = RegisterObject("Caribou hyperelastic force field")
.add<HyperelasticForcefieldStoreF<Hexahedron<Linear>>>()
.add<HyperelasticForcefieldStoreF<Tetrahedron<Linear>>>();

static int _c3_ = RegisterObject("Caribou hyperelastic force field")
.add<HyperelasticForcefieldStoreFAndS<Hexahedron<Linear>>>()
.add<HyperelasticForcefieldStoreFAndS<Tetrahedron<Linear>>>();