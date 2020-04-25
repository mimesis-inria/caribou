#include <sofa/core/ObjectFactory.h>
#include <SofaCaribou/Topology/CircleIsoSurface.h>
#include <SofaCaribou/Topology/SphereIsoSurface.h>
#include <SofaCaribou/Topology/CylinderIsoSurface.h>

namespace SofaCaribou::topology {
using namespace sofa::core;

static int CircleIsoSurfaceClass = RegisterObject("Caribou circle iso-surface.").add<CircleIsoSurface>(true);
static int SphereIsoSurfaceClass = RegisterObject("Caribou sphere iso-surface.").add<SphereIsoSurface>(true);
static int CylinderIsoSurfaceClass = RegisterObject("Caribou cylinder iso-surface.").add<CylinderIsoSurface>(true);

}
