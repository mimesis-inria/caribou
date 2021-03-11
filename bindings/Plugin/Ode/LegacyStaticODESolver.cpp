#include "StaticODESolver.h"

#include <SofaCaribou/Ode/LegacyStaticODESolver.h>

#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include <SofaPython3/PythonFactory.h>
#include <SofaPython3/Sofa/Core/Binding_Base.h>

namespace py = pybind11;

PYBIND11_MAKE_OPAQUE(std::vector<FLOATING_POINT_TYPE>)
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<FLOATING_POINT_TYPE>>)

namespace SofaCaribou::ode::python {

void addLegacyStaticODESolver(py::module &m) {
    using namespace sofa::core::objectmodel;
    py::class_<LegacyStaticODESolver, sofa::core::objectmodel::BaseObject, sofapython3::py_shared_ptr<LegacyStaticODESolver>> c (m, "LegacyStaticODESolver");
    c.def_property_readonly("iteration_times", &LegacyStaticODESolver::iteration_times);
    c.def_property_readonly("squared_residuals", &LegacyStaticODESolver::squared_residuals);
    c.def_property_readonly("squared_initial_residual", &LegacyStaticODESolver::squared_initial_residual);
    c.def_property_readonly("iterative_linear_solver_squared_residuals", &LegacyStaticODESolver::iterative_linear_solver_squared_residuals);
    c.def_property_readonly("iterative_linear_solver_squared_rhs_norms", &LegacyStaticODESolver::iterative_linear_solver_squared_rhs_norms);

    py::bind_vector<std::vector<FLOATING_POINT_TYPE>>(m, "VectorFloat");
    py::bind_vector<std::vector<std::vector<FLOATING_POINT_TYPE>>>(m, "VectorVectorFloat");

    sofapython3::PythonFactory::registerType<LegacyStaticODESolver>([](sofa::core::objectmodel::Base* o) {
        return py::cast(dynamic_cast<LegacyStaticODESolver*>(o));
    });
}

} // namespace SofaCaribou::ode::python