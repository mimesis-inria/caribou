#pragma once

#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace SofaCaribou::solver::python {

void addConjugateGradientSolver(py::module &m);

}