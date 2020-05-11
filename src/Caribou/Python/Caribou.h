#pragma once

#include <pybind11/pybind11.h>
#include <Caribou/Constants.h>

namespace caribou::python {

enum Dimension {
    _1D = caribou::_1D,
    _2D = caribou::_2D,
    _3D = caribou::_3D
};

enum Order {
    Linear = caribou::Linear,
    Quadratic = caribou::Quadratic
};

}