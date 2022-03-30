#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>

#include <SofaPython3/PythonFactory.h>
#include <SofaPython3/Sofa/Core/Binding_Base.h>

#include <SofaCaribou/FEniCS/Forcefield/HyperelasticForcefield_FEniCS.h>
#include <SofaCaribou/FEniCS/Forcefield/HyperelasticForcefield_FEniCS.inl>

namespace SofaCaribou::forcefield::python {

template<typename Element>
void bind_hyperelastic_forcefield_FEniCS(pybind11::module &m, const std::string & template_name) {
    pybind11::module::import("Sofa");

    std::string name = "HyperelasticForcefield_FEniCS<" + template_name + ">";

    pybind11::class_<HyperelasticForcefield_FEniCS<Element>, sofa::core::objectmodel::BaseObject, sofapython3::py_shared_ptr<HyperelasticForcefield_FEniCS<Element>>> c(m, name.c_str());

    c.def("K", &HyperelasticForcefield_FEniCS<Element>::K);
    c.def("energy", &HyperelasticForcefield_FEniCS<Element>::getPotentialEnergy);
    c.def("cond", &HyperelasticForcefield_FEniCS<Element>::cond);
    c.def("eigenvalues", &HyperelasticForcefield_FEniCS<Element>::eigenvalues);
    c.def("assemble_stiffness", [](HyperelasticForcefield_FEniCS<Element> & self, const Eigen::Matrix<double, Eigen::Dynamic, HyperelasticForcefield_FEniCS<Element>::Dimension, Eigen::RowMajor> & x, const Eigen::Matrix<double, Eigen::Dynamic, HyperelasticForcefield_FEniCS<Element>::Dimension, Eigen::RowMajor> & x0) {
        self.assemble_stiffness(x, x0);
    }, pybind11::arg("x").noconvert(true), pybind11::arg("x0").noconvert(true));
    c.def("assemble_stiffness", [](HyperelasticForcefield_FEniCS<Element> & self, const Eigen::Matrix<float, Eigen::Dynamic, HyperelasticForcefield_FEniCS<Element>::Dimension, Eigen::RowMajor> & x, const Eigen::Matrix<float, Eigen::Dynamic, HyperelasticForcefield_FEniCS<Element>::Dimension, Eigen::RowMajor> & x0) {
        self.assemble_stiffness(x, x0);
    }, pybind11::arg("x").noconvert(true), pybind11::arg("x").noconvert(true));
    sofapython3::PythonFactory::registerType<HyperelasticForcefield_FEniCS<Element>>([template_name](sofa::core::objectmodel::Base* o) {
        return pybind11::cast(dynamic_cast<HyperelasticForcefield_FEniCS<Element>*>(o));
    });
}

template<typename Element>
void bind_hyperelastic_forcefield_FEniCS(pybind11::module &m) {
    bind_hyperelastic_forcefield_FEniCS<Element>(m, HyperelasticForcefield_FEniCS<Element>::templateName());
}

void addHyperElasticForcefield_FEniCS(pybind11::module &m);
}
