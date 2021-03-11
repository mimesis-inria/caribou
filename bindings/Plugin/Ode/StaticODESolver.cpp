#include "StaticODESolver.h"

#include <SofaCaribou/Ode/StaticODESolver.h>

#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include <SofaPython3/PythonFactory.h>
#include <SofaPython3/Sofa/Core/Binding_Base.h>

namespace py = pybind11;

PYBIND11_MAKE_OPAQUE(std::vector<FLOATING_POINT_TYPE>)
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<FLOATING_POINT_TYPE>>)

namespace SofaCaribou::ode::python {

void addStaticODESolver(py::module &m) {
    using namespace sofa::core::objectmodel;
    py::class_<StaticODESolver, sofa::core::objectmodel::BaseObject, sofapython3::py_shared_ptr<StaticODESolver>> c (m, "StaticODESolver");
    c.def_property_readonly("iteration_times", &StaticODESolver::iteration_times);
    c.def_property_readonly("squared_residuals", &StaticODESolver::squared_residuals);
    c.def_property_readonly("squared_initial_residual", &StaticODESolver::squared_initial_residual);

    py::bind_vector<std::vector<FLOATING_POINT_TYPE>>(m, "VectorFloat");
    py::bind_vector<std::vector<std::vector<FLOATING_POINT_TYPE>>>(m, "VectorVectorFloat");

    sofapython3::PythonFactory::registerType<StaticODESolver>([](sofa::core::objectmodel::Base* o) {
        return py::cast(dynamic_cast<StaticODESolver*>(o));
    });
}

} // namespace SofaCaribou::ode::python