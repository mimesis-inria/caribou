#include "ConjugateGradientSolver.h"

#include <SofaCaribou/Solver/ConjugateGradientSolver.h>

#include <pybind11/eigen.h>

#include <SofaPython3/PythonFactory.h>
#include <SofaPython3/Sofa/Core/Binding_Base.h>

namespace py = pybind11;

namespace SofaCaribou::solver::python {

void addConjugateGradientSolver(py::module & m) {
    py::class_<ConjugateGradientSolver, sofa::core::objectmodel::BaseObject, sofapython3::py_shared_ptr<ConjugateGradientSolver>> c(m, "ConjugateGradientSolver");

    c.def("K", &ConjugateGradientSolver::system_matrix);

    c.def("assemble", [](ConjugateGradientSolver & solver, double m, double b, double k) {
        sofa::core::MechanicalParams mparams;
        mparams.setMFactor(m);
        mparams.setBFactor(b);
        mparams.setKFactor(k);
        solver.assemble(&mparams);
    }, py::arg("m") = static_cast<double>(1), py::arg("b") = static_cast<double>(1), py::arg("k") = static_cast<double>(1));

    sofapython3::PythonFactory::registerType<ConjugateGradientSolver>([](sofa::core::objectmodel::Base* o) {
        return py::cast(dynamic_cast<ConjugateGradientSolver*>(o));
    });
}

} // namespace SofaCaribou::solver::python