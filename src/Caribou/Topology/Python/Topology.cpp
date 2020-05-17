#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace caribou::topology::python {
void create_hashgrid(pybind11::module & m);
void create_grid(pybind11::module & m);
void create_mesh(py::module & m);
}

namespace caribou::topology::io::python {
void create_IO(pybind11::module &m);
}

PYBIND11_MODULE(CaribouTopologyPython, m) {
    m.doc() = "Topology module";
    caribou::topology::python::create_grid(m);
    caribou::topology::python::create_hashgrid(m);
    caribou::topology::python::create_mesh(m);
    caribou::topology::io::python::create_IO(m);
}