#pragma once

namespace pybind11{ class module; }

namespace caribou::geometry::python {

void create_segment(pybind11::module & m);

}