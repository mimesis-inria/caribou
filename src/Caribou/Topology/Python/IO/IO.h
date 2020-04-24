#pragma once

namespace pybind11 {class module;}

namespace caribou::topology::io::python {

void create_IO(pybind11::module & m);

} // namespace caribou::topology::io::python