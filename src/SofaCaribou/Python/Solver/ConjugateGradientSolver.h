#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>

#include <SofaCaribou/Solver/ConjugateGradientSolver.h>

#include <SofaPython3/PythonFactory.h>
#include <SofaPython3/Sofa/Core/Binding_Base.h>

namespace SofaCaribou::solver::python {

template <typename EigenMatrix>
void bind_ConjugateGradientSolver(pybind11::module & m) {
    namespace py = pybind11;
    using SolverType = ConjugateGradientSolver<EigenMatrix>;
    py::class_<SolverType, sofa::core::objectmodel::BaseObject, sofapython3::py_shared_ptr<SolverType>> c(m, "ConjugateGradientSolver");

    c.def("A", [](const SolverType & solver){return solver.A()->matrix();});

    c.def("x", [](const SolverType & solver){return solver.x()->vector();});

    c.def("b", [](const SolverType & solver){return solver.b()->vector();});


    c.def("assemble", [](SolverType & solver, double m, double b, double k) {
        sofa::core::MechanicalParams mparams;
        mparams.setMFactor(m);
        mparams.setBFactor(b);
        mparams.setKFactor(k);
        solver.assemble(&mparams);
    }, py::arg("m") = static_cast<double>(1), py::arg("b") = static_cast<double>(1), py::arg("k") = static_cast<double>(1));

    sofapython3::PythonFactory::registerType<SolverType>([](sofa::core::objectmodel::Base* o) {
        return py::cast(dynamic_cast<SolverType*>(o));
    });
}

void addConjugateGradientSolver(pybind11::module &m);

}