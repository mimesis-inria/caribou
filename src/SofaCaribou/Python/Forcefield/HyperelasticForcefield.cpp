#include "HyperelasticForcefield.h"

#include <pybind11/eigen.h>

namespace SofaCaribou::forcefield::python {

void addHyperElasticForcefield(py::module &m) {
    using namespace caribou;
    using namespace caribou::geometry;
    bind_hyperelastic_forcefield<Tetrahedron<Linear>>(m);
    bind_hyperelastic_forcefield<Hexahedron<Linear>>(m);
}
}