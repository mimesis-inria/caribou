#include <pybind11/pybind11.h>

#include <SofaCaribou/Python/Topology/FictitiousGrid.h>
#include <SofaCaribou/Python/Forcefield/FictitiousGridElasticForce.h>
#include <SofaCaribou/Python/Forcefield/HexahedronElasticForce.h>

namespace py = pybind11;


PYBIND11_MODULE(SofaCaribouPython, m) {
    m.doc() = "SofaCaribou module";

    SofaCaribou::topology::python::addFictitiousGrid(m);
    SofaCaribou::forcefield::python::addHexahedronElasticForce(m);
    SofaCaribou::forcefield::python::addFictitiousGridElasticForce(m);
}