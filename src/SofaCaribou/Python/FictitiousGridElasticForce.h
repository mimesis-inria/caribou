#pragma once

#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace SofaCaribou::Python {
void addFictitiousGridElasticForce(py::module &m);
}