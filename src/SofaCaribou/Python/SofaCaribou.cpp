#include <pybind11/pybind11.h>

#include <SofaCaribou/Python/Topology/FictitiousGrid.h>
#include <SofaCaribou/Python/Forcefield/FictitiousGridElasticForce.h>
#include <SofaCaribou/Python/Ode/StaticODESolver.h>
#include <SofaCaribou/Python/Forcefield/HexahedronElasticForce.h>
#include <SofaCaribou/Python/Forcefield/HyperelasticForcefield.h>
#include <SofaCaribou/Python/Forcefield/FictitiousGridHyperelasticForcefield.h>

namespace py = pybind11;


PYBIND11_MODULE(SofaCaribouPython, m) {
    m.doc() = "SofaCaribou module";

    // Topology bindings
    SofaCaribou::topology::python::addFictitiousGrid(m);

    // ODE bindings
    SofaCaribou::ode::python::addStaticODESolver(m);

    // Forcefield bindings
    SofaCaribou::forcefield::python::addFictitiousGridElasticForce(m);
    SofaCaribou::forcefield::python::addFictitiousGridHyperelasticForcefield(m);
    SofaCaribou::forcefield::python::addHexahedronElasticForce(m);
    SofaCaribou::forcefield::python::addHyperElasticForcefield(m);
}