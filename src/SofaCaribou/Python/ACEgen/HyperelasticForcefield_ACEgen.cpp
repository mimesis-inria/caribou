#include "HyperelasticForcefield_ACEgen.h"

#include <Caribou/constants.h>
#include <Caribou/Geometry/Tetrahedron.h>
#include <Caribou/Geometry/Hexahedron.h>

#include <pybind11/eigen.h>

namespace SofaCaribou::forcefield::python {

void addHyperElasticForcefield_ACEgen(pybind11::module &m) {
    using namespace caribou;
    using namespace caribou::geometry;
//    bind_hyperelastic_forcefield_ACEgen<Tetrahedron<Linear>>(m);
//    bind_hyperelastic_forcefield_ACEgen<Tetrahedron<Quadratic>>(m);
    bind_hyperelastic_forcefield_ACEgen<Hexahedron<Linear>>(m);
//    bind_hyperelastic_forcefield_ACEgen<Hexahedron<Quadratic>>(m);
}

}
