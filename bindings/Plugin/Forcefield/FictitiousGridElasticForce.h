#pragma once

#include <pybind11/pybind11.h>

namespace SofaCaribou::forcefield::python {
void addFictitiousGridElasticForce(pybind11::module &m);
}