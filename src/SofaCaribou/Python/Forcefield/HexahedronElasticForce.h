#pragma once

#include <pybind11/pybind11.h>

namespace SofaCaribou::forcefield::python {
void addHexahedronElasticForce(pybind11::module &m);
}