#include "HexahedronElasticForce.h"
#include <SofaCaribou/GraphComponents/Forcefield/HexahedronElasticForce.h>

namespace SofaCaribou::Python {
void addHexahedronElasticForce(py::module &m) {
    //    py::object ForceField = (py::object) py::module::import("Sofa.Core").attr("ForceField");
    using HexahedronElasticForce = SofaCaribou::GraphComponents::forcefield::HexahedronElasticForce;
    py::class_<HexahedronElasticForce> c(m, "HexahedronElasticForce");
}
}