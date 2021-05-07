#include <pybind11/pybind11.h>

#include <SofaCaribou/Python/Ode/LegacyStaticODESolver.h>
#include <SofaCaribou/Python/Ode/StaticODESolver.h>
#include <SofaCaribou/Python/Forcefield/HexahedronElasticForce.h>
#include <SofaCaribou/Python/Forcefield/HyperelasticForcefield.h>
#include <SofaCaribou/Python/Solver/ConjugateGradientSolver.h>
#include <SofaCaribou/Python/Topology/CaribouTopology.h>
#include <SofaCaribou/Python/Topology/FictitiousGrid.h>

#include <vector>
#include <pybind11/stl_bind.h>

PYBIND11_MODULE(SofaCaribou, m) {
    m.doc() = "SofaCaribou module";

    // Topology bindings
    SofaCaribou::topology::python::addFictitiousGrid(m);
    SofaCaribou::topology::python::addCaribouTopology(m);

    // ODE bindings
    SofaCaribou::ode::python::addLegacyStaticODESolver(m);
    SofaCaribou::ode::python::addStaticODESolver(m);

    // Forcefield bindings
    SofaCaribou::forcefield::python::addHexahedronElasticForce(m);
    SofaCaribou::forcefield::python::addHyperElasticForcefield(m);

    // Solver bindings
    SofaCaribou::solver::python::addConjugateGradientSolver(m);

    // Container bindings
    pybind11::bind_vector<std::vector<FLOATING_POINT_TYPE>>(m, "VectorFloat");
    pybind11::bind_vector<std::vector<std::vector<FLOATING_POINT_TYPE>>>(m, "VectorVectorFloat");
}