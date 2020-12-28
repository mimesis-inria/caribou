#pragma once

#include <pybind11/pybind11.h>

namespace SofaCaribou::forcefield::python {
void addFictitiousGridHyperelasticForcefield(pybind11::module &m);
}