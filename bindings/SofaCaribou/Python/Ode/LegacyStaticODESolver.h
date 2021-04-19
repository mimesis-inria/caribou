#pragma once

#include <pybind11/pybind11.h>

namespace SofaCaribou::ode::python {

void addLegacyStaticODESolver(pybind11::module &m);

}