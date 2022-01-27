#include "HyperelasticForcefield_FEniCS.h"

#include <Caribou/constants.h>
#include <Caribou/Geometry/Tetrahedron.h>
#include <Caribou/Geometry/Hexahedron.h>

#include <pybind11/eigen.h>

namespace SofaCaribou::forcefield::python {

void addHyperElasticForcefield_FEniCS(pybind11::module &m) {
    using namespace caribou;
    using namespace caribou::geometry;
    bind_hyperelastic_forcefield_FEniCS<Tetrahedron<Linear>>(m);
    bind_hyperelastic_forcefield_FEniCS<Tetrahedron<Quadratic>>(m);
    bind_hyperelastic_forcefield_FEniCS<Hexahedron<Linear>>(m);
    bind_hyperelastic_forcefield_FEniCS<Hexahedron<Quadratic>>(m);
}

}
