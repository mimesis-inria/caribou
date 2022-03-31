#include "HyperelasticForcefield_FEniCS.h"

#include <Caribou/constants.h>
#include <Caribou/Geometry/Tetrahedron.h>
#include <Caribou/Geometry/Tetrahedron10.h>
#include <Caribou/Geometry/Hexahedron_FEniCS.h>
#include <Caribou/Geometry/Hexahedron_FEniCS20.h>

#include <pybind11/eigen.h>

namespace SofaCaribou::forcefield::python {

void addHyperElasticForcefield_FEniCS(pybind11::module &m) {
    using namespace caribou;
    using namespace caribou::geometry;
    bind_hyperelastic_forcefield_FEniCS<Tetrahedron>(m);
    bind_hyperelastic_forcefield_FEniCS<Tetrahedron10>(m);
    bind_hyperelastic_forcefield_FEniCS<Hexahedron_FEniCS>(m);
    bind_hyperelastic_forcefield_FEniCS<Hexahedron_FEniCS20>(m);

}

}
