#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>

#include <SofaCaribou/Solver/LDLTSolver.h>
#include<Eigen/SparseCholesky>

#include <SofaPython3/PythonFactory.h>
#include <SofaPython3/Sofa/Core/Binding_Base.h>

namespace SofaCaribou::solver::python {

    template <typename EigenSolver>
    void bind_LDLTSolver(pybind11::module & m) {
        namespace py = pybind11;
        using SOLVER = SofaCaribou::solver::LDLTSolver<EigenSolver>;
        py::class_<SOLVER, sofa::core::objectmodel::BaseObject, sofapython3::py_shared_ptr<SOLVER>> c(m, "LDLTSolver");

        c.def("A", [](SOLVER & solver){return solver.A()->matrix();});

        c.def("x", [](SOLVER & solver){return solver.x()->vector();});

        c.def("b", [](SOLVER & solver){return solver.b()->vector();});


        c.def("assemble", [](SOLVER & solver, double m, double b, double k) {
            sofa::core::MechanicalParams mparams;
            mparams.setMFactor(m);
            mparams.setBFactor(b);
            mparams.setKFactor(k);
            solver.assemble(&mparams);
        }, py::arg("m") = static_cast<double>(1), py::arg("b") = static_cast<double>(1), py::arg("k") = static_cast<double>(1));

        sofapython3::PythonFactory::registerType<SOLVER>([](sofa::core::objectmodel::Base* o) {
            return py::cast(dynamic_cast<SOLVER*>(o));
        });
    }

    void addLDLTSolver(pybind11::module &m);

}