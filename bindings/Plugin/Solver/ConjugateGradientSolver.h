#pragma once

#include <pybind11/pybind11.h>

namespace SofaCaribou::solver::python {

void addConjugateGradientSolver(pybind11::module &m);

}