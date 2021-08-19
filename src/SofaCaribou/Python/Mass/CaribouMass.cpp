#include "CaribouMass.h"

#include <Caribou/constants.h>
#include <Caribou/Geometry/Tetrahedron.h>
#include <Caribou/Geometry/Hexahedron.h>
#include <SofaCaribou/Mass/CaribouMass[Tetrahedron].h>
#include <SofaCaribou/Mass/CaribouMass[Hexahedron].h>
#include <pybind11/eigen.h>

namespace SofaCaribou::mass::python {

void addCaribouMass(pybind11::module &m) {
    using namespace caribou;
    using namespace caribou::geometry;
    bind_caribou_mass<Tetrahedron<Linear>>(m);
    bind_caribou_mass<Tetrahedron<Quadratic>>(m);
    bind_caribou_mass<Hexahedron<Linear>>(m);
    bind_caribou_mass<Hexahedron<Quadratic>>(m);
}
}