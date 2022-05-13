#include "HyperelasticForcefield.h"

#include <Caribou/constants.h>
#include <Caribou/Geometry/Tetrahedron.h>
#include <Caribou/Geometry/Tetrahedron10.h>
#include <Caribou/Geometry/Hexahedron.h>
#include <Caribou/Geometry/Hexahedron20.h>

#include <pybind11/eigen.h>

namespace SofaCaribou::forcefield::python {

void addHyperElasticForcefield(pybind11::module &m) {
    using namespace caribou;
    using namespace caribou::geometry;
    bind_hyperelastic_forcefield<Tetrahedron>(m);
    bind_hyperelastic_forcefield<Tetrahedron10>(m);
    bind_hyperelastic_forcefield<Hexahedron>(m);
    bind_hyperelastic_forcefield<Hexahedron20>(m);
}

}