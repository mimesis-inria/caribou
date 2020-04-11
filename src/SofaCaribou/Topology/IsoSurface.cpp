#include <SofaCaribou/Topology/CircleIsoSurface.h>
#include <SofaCaribou/Topology/SphereIsoSurface.h>
#include <sofa/core/ObjectFactory.h>

namespace SofaCaribou::topology {
using namespace sofa::core;

static int CircleIsoSurfaceClass = RegisterObject("Caribou circle iso-surface.").add<CircleIsoSurface>(true);
static int SphereIsoSurfaceClass = RegisterObject("Caribou sphere iso-surface.").add<SphereIsoSurface>(true);

}