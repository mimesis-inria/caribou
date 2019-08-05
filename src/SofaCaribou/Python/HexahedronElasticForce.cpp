#include "HexahedronElasticForce.h"
#include <SofaCaribou/GraphComponents/Forcefield/HexahedronElasticForce.h>

#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>

PYBIND11_MAKE_OPAQUE(const std::vector<SofaCaribou::GraphComponents::forcefield::HexahedronElasticForce::GaussNode>&);

namespace SofaCaribou::Python {
void addHexahedronElasticForce(py::module &m) {
    //    py::object ForceField = (py::object) py::module::import("Sofa.Core").attr("ForceField");
    using HexahedronElasticForce = SofaCaribou::GraphComponents::forcefield::HexahedronElasticForce;

    py::class_<HexahedronElasticForce::GaussNode> g(m, "GaussNode");
    g.def_property_readonly("weight", [](const HexahedronElasticForce::GaussNode & self){return self.weight;});
    g.def_property_readonly("jacobian_determinant", [](const HexahedronElasticForce::GaussNode & self){return self.jacobian_determinant;});
    g.def_property_readonly("dN_dx", [](const HexahedronElasticForce::GaussNode & self){return self.dN_dx;});
    g.def_property_readonly("F", [](const HexahedronElasticForce::GaussNode & self){return self.F;});

    py::class_<HexahedronElasticForce> c(m, "HexahedronElasticForce");
    c.def("gauss_nodes_of", &HexahedronElasticForce::gauss_nodes_of, py::arg("hexahedron_id"), py::return_value_policy::reference_internal);
    c.def("stiffness_matrix_of", &HexahedronElasticForce::stiffness_matrix_of, py::arg("hexahedron_id"), py::return_value_policy::reference_internal);
}
}