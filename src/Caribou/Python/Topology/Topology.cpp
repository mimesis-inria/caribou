#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace caribou::topology::bindings {
void create_hashgrid(pybind11::module & m);
void create_grid(pybind11::module & m);
void create_mesh(py::module & m);
void create_domain(py::module & m);
}

namespace caribou::topology::io::bindings {
void create_IO(pybind11::module &m);
}

PYBIND11_MODULE(Topology, m) {
    m.doc() = "Topology module";
    caribou::topology::bindings::create_grid(m);
    caribou::topology::bindings::create_hashgrid(m);
    caribou::topology::bindings::create_domain(m);
    caribou::topology::bindings::create_mesh(m);
    caribou::topology::io::bindings::create_IO(m);
}