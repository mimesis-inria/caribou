#include "HexahedronElasticForce.h"
#include <SofaCaribou/Forcefield/FictitiousGridElasticForce.h>

#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <memory>

PYBIND11_MAKE_OPAQUE(const std::vector<SofaCaribou::GraphComponents::forcefield::FictitiousGridElasticForce::GaussNode>&)

using sofa::defaulttype::Vec3Types;
using namespace sofa::core::objectmodel;
using namespace sofa::core::behavior;

namespace SofaCaribou::Python {
void addFictitiousGridElasticForce(py::module &m) {
    using FictitiousGridElasticForce = SofaCaribou::GraphComponents::forcefield::FictitiousGridElasticForce;

    py::class_<FictitiousGridElasticForce::GaussNode> g(m, "FictitiousGridElasticForceGaussNode");
    g.def_property_readonly("weight", [](const FictitiousGridElasticForce::GaussNode & self){return self.weight;});
    g.def_property_readonly("dN_dx", [](const FictitiousGridElasticForce::GaussNode & self){return self.dN_dx;});
    g.def_property_readonly("F", [](const FictitiousGridElasticForce::GaussNode & self){return self.F;});

//    py::object BaseObject = (py::object) py::module::import("Sofa.Core").attr("Base");
    py::class_<FictitiousGridElasticForce> c(m, "FictitiousGridElasticForce");
    c.def("gauss_nodes_of", &FictitiousGridElasticForce::gauss_nodes_of, py::arg("hexahedron_id"), py::return_value_policy::reference_internal);
    c.def("stiffness_matrix_of", &FictitiousGridElasticForce::stiffness_matrix_of, py::arg("hexahedron_id"), py::return_value_policy::reference_internal);
    c.def("K", &FictitiousGridElasticForce::K);
    c.def("cond", &FictitiousGridElasticForce::cond);
    c.def("eigenvalues", &FictitiousGridElasticForce::eigenvalues);
}
}