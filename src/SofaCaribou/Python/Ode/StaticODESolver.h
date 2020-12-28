#pragma once

#include <pybind11/pybind11.h>

namespace SofaCaribou::ode::python {

void addStaticODESolver(pybind11::module &m);

}