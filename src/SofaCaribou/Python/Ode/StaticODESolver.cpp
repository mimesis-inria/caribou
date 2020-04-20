#include "StaticODESolver.h"

#include <SofaCaribou/Ode/StaticODESolver.h>

#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

PYBIND11_MAKE_OPAQUE(std::vector<FLOATING_POINT_TYPE>)
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<FLOATING_POINT_TYPE>>)

namespace SofaCaribou::ode::python {

void addStaticODESolver(py::module &m) {
    using namespace sofa::core::objectmodel;
    py::class_<StaticODESolver, std::shared_ptr<StaticODESolver>> c (m, "StaticODESolver");
    c.def_property_readonly("iteration_times", &StaticODESolver::iteration_times);
    c.def_property_readonly("squared_residuals", &StaticODESolver::squared_residuals);
    c.def_property_readonly("squared_initial_residual", &StaticODESolver::squared_initial_residual);
    c.def_property_readonly("iterative_linear_solver_squared_residuals", &StaticODESolver::iterative_linear_solver_squared_residuals);
    c.def_property_readonly("iterative_linear_solver_squared_rhs_norms", &StaticODESolver::iterative_linear_solver_squared_rhs_norms);

    py::bind_vector<std::vector<FLOATING_POINT_TYPE>>(m, "VectorFloat");
    py::bind_vector<std::vector<std::vector<FLOATING_POINT_TYPE>>>(m, "VectorVectorFloat");
}

} // namespace SofaCaribou::ode::python