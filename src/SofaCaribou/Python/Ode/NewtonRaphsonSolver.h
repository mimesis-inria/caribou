#pragma once

#undef NDEBUG
#include <pybind11/pybind11.h>

namespace SofaCaribou::ode::python {

void addNewtonRaphsonSolver(pybind11::module &m);

}
