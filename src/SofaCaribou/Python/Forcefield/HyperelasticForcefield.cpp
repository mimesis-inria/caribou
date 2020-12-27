#include "HyperelasticForcefield.h"

#include <Caribou/Constants.h>
#include <Caribou/Geometry/Tetrahedron.h>
#include <Caribou/Geometry/Hexahedron.h>

#include <pybind11/eigen.h>

namespace SofaCaribou::forcefield::python {

void addHyperElasticForcefield(pybind11::module &m) {
    using namespace caribou;
    using namespace caribou::geometry;
    bind_hyperelastic_forcefield<Tetrahedron<Linear>>(m);
    bind_hyperelastic_forcefield<Hexahedron<Linear>>(m);
}
}