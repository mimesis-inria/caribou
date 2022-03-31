#include "CaribouMass.h"

#include <Caribou/constants.h>
#include <Caribou/Geometry/Tetrahedron.h>
#include <Caribou/Geometry/Tetrahedron10.h>
#include <Caribou/Geometry/Hexahedron.h>
#include <Caribou/Geometry/Hexahedron20.h>
#include <SofaCaribou/Mass/CaribouMass[Tetrahedron].h>
#include <SofaCaribou/Mass/CaribouMass[Tetrahedron10].h>
#include <SofaCaribou/Mass/CaribouMass[Hexahedron].h>
#include <SofaCaribou/Mass/CaribouMass[Hexahedron20].h>
#include <pybind11/eigen.h>

namespace SofaCaribou::mass::python {

void addCaribouMass(pybind11::module &m) {
    using namespace caribou;
    using namespace caribou::geometry;
    bind_caribou_mass<Tetrahedron>(m);
    bind_caribou_mass<Tetrahedron10>(m);
    bind_caribou_mass<Hexahedron>(m);
    bind_caribou_mass<Hexahedron20>(m);
}
}