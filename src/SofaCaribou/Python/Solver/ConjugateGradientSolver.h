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
    py::class_<ConjugateGradientSolver<EigenMatrix>, sofa::core::objectmodel::BaseObject, sofapython3::py_shared_ptr<ConjugateGradientSolver<EigenMatrix>>> c(m, "ConjugateGradientSolver");

    c.def("A", &ConjugateGradientSolver<EigenMatrix>::A);

    c.def("assemble", [](ConjugateGradientSolver<EigenMatrix> & solver, double m, double b, double k) {
        sofa::core::MechanicalParams mparams;
        mparams.setMFactor(m);
        mparams.setBFactor(b);
        mparams.setKFactor(k);
        solver.assemble(&mparams);
    }, py::arg("m") = static_cast<double>(1), py::arg("b") = static_cast<double>(1), py::arg("k") = static_cast<double>(1));

    sofapython3::PythonFactory::registerType<ConjugateGradientSolver<EigenMatrix>>([](sofa::core::objectmodel::Base* o) {
        return py::cast(dynamic_cast<ConjugateGradientSolver<EigenMatrix>*>(o));
    });
}

void addConjugateGradientSolver(pybind11::module &m);

}