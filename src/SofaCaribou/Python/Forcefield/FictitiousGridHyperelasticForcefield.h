#pragma once

#include <pybind11/pybind11.h>

#include <SofaCaribou/Python/Forcefield/HyperelasticForcefield.h>

namespace SofaCaribou::forcefield::python {
void addFictitiousGridHyperelasticForcefield(py::module &m);
}