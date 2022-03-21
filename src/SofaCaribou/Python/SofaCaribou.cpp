#include <pybind11/pybind11.h>

#include <SofaCaribou/Python/Ode/LegacyStaticODESolver.h>
#include <SofaCaribou/Python/Ode/StaticODESolver.h>
#include <SofaCaribou/Python/Mass/CaribouMass.h>
#include <SofaCaribou/Python/Forcefield/HexahedronElasticForce.h>
#include <SofaCaribou/Python/Forcefield/HyperelasticForcefield.h>
#include <SofaCaribou/Python/FEniCS/HyperelasticForcefield_FEniCS.h>
#include <SofaCaribou/Python/ACEgen/HyperelasticForcefield_ACEgen.h>
#include <SofaCaribou/Python/Solver/ConjugateGradientSolver.h>
#include <SofaCaribou/Python/Solver/LDLTSolver.h>
#include <SofaCaribou/Python/Solver/LLTSolver.h>
#include <SofaCaribou/Python/Solver/LUSolver.h>
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

    // Mass bindings
    SofaCaribou::mass::python::addCaribouMass(m);

    // Forcefield bindings
    SofaCaribou::forcefield::python::addHexahedronElasticForce(m);
    SofaCaribou::forcefield::python::addHyperElasticForcefield(m);
    SofaCaribou::forcefield::python::addHyperElasticForcefield_FEniCS(m);
//    SofaCaribou::forcefield::python::addHyperElasticForcefield_ACEgen(m);
    // Solver bindings
    SofaCaribou::solver::python::addConjugateGradientSolver(m);
    SofaCaribou::solver::python::addLDLTSolver(m);
    SofaCaribou::solver::python::addLLTSolver(m);
    SofaCaribou::solver::python::addLUSolver(m);


    // Container bindings
    pybind11::bind_vector<std::vector<FLOATING_POINT_TYPE>>(m, "VectorFloat");
    pybind11::bind_vector<std::vector<std::vector<FLOATING_POINT_TYPE>>>(m, "VectorVectorFloat");
}
