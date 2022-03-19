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
    py::class_<StaticODESolver, NewtonRaphsonSolver, sofapython3::py_shared_ptr<StaticODESolver>> c (m, "StaticODESolver");

    sofapython3::PythonFactory::registerType<StaticODESolver>([](sofa::core::objectmodel::Base* o) {
        return py::cast(dynamic_cast<StaticODESolver*>(o));
    });
}

} // namespace SofaCaribou::ode::python