#include "HexahedronElasticForce.h"
#include <SofaCaribou/Forcefield/HexahedronElasticForce.h>
#include <SofaPython3/Sofa/Core/Binding_Base.h>
#include <SofaPython3/PythonFactory.h>

#include <pybind11/stl.h>
#include <pybind11/eigen.h>

PYBIND11_MAKE_OPAQUE(const std::vector<SofaCaribou::forcefield::HexahedronElasticForce::GaussNode>&)

using sofa::defaulttype::Vec3Types;
using namespace sofa::core::objectmodel;
using namespace sofa::core::behavior;
namespace py = pybind11;


namespace SofaCaribou::forcefield::python {

void addHexahedronElasticForce(py::module &m) {
    using HexahedronElasticForce = SofaCaribou::forcefield::HexahedronElasticForce;

    py::class_<HexahedronElasticForce::GaussNode> g(m, "HexahedronElasticForceGaussNode");
    g.def_property_readonly("weight", [](const HexahedronElasticForce::GaussNode & self){return self.weight;});
    g.def_property_readonly("jacobian_determinant", [](const HexahedronElasticForce::GaussNode & self){return self.jacobian_determinant;});
    g.def_property_readonly("dN_dx", [](const HexahedronElasticForce::GaussNode & self){return self.dN_dx;});
    g.def_property_readonly("F", [](const HexahedronElasticForce::GaussNode & self){return self.F;});

    py::class_<HexahedronElasticForce, sofa::core::objectmodel::BaseObject, sofapython3::py_shared_ptr<HexahedronElasticForce>> c(m, "HexahedronElasticForce");
    c.def("gauss_nodes_of", &HexahedronElasticForce::gauss_nodes_of, py::arg("hexahedron_id"), py::return_value_policy::reference_internal);
    c.def("stiffness_matrix_of", &HexahedronElasticForce::stiffness_matrix_of, py::arg("hexahedron_id"), py::return_value_policy::reference_internal);
    c.def("K", &HexahedronElasticForce::K);
    c.def("cond", &HexahedronElasticForce::cond);
    c.def("eigenvalues", &HexahedronElasticForce::eigenvalues);

    sofapython3::PythonFactory::registerType<HexahedronElasticForce>([](sofa::core::objectmodel::Base* o) {
        return pybind11::cast(dynamic_cast<HexahedronElasticForce*>(o));
    });
}

}