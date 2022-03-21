#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>

#include <SofaPython3/PythonFactory.h>
#include <SofaPython3/Sofa/Core/Binding_Base.h>

#include <SofaCaribou/ACEgen/Forcefield/HyperelasticForcefield_ACEgen.h>
#include <SofaCaribou/ACEgen/Forcefield/HyperelasticForcefield_ACEgen.inl>

namespace SofaCaribou::forcefield::python {

template<typename Element>
void bind_hyperelastic_forcefield_ACEgen(pybind11::module &m, const std::string & template_name) {
    pybind11::module::import("Sofa");

    std::string name = "HyperelasticForcefield_ACEgen<" + template_name + ">";

    pybind11::class_<HyperelasticForcefield_ACEgen<Element>, sofa::core::objectmodel::BaseObject, sofapython3::py_shared_ptr<HyperelasticForcefield_ACEgen<Element>>> c(m, name.c_str());

    c.def("K", &HyperelasticForcefield_ACEgen<Element>::K);
    c.def("cond", &HyperelasticForcefield_ACEgen<Element>::cond);
    c.def("eigenvalues", &HyperelasticForcefield_ACEgen<Element>::eigenvalues);
    c.def("assemble_stiffness", [](HyperelasticForcefield_ACEgen<Element> & self, const Eigen::Matrix<double, Eigen::Dynamic, HyperelasticForcefield_ACEgen<Element>::Dimension, Eigen::RowMajor> & x, const Eigen::Matrix<double, Eigen::Dynamic, HyperelasticForcefield_ACEgen<Element>::Dimension, Eigen::RowMajor> & x0) {
        self.assemble_stiffness(x, x0);
    }, pybind11::arg("x").noconvert(true), pybind11::arg("x0").noconvert(true));
    c.def("assemble_stiffness", [](HyperelasticForcefield_ACEgen<Element> & self, const Eigen::Matrix<float, Eigen::Dynamic, HyperelasticForcefield_ACEgen<Element>::Dimension, Eigen::RowMajor> & x, const Eigen::Matrix<float, Eigen::Dynamic, HyperelasticForcefield_ACEgen<Element>::Dimension, Eigen::RowMajor> & x0) {
        self.assemble_stiffness(x, x0);
    }, pybind11::arg("x").noconvert(true), pybind11::arg("x").noconvert(true));
    sofapython3::PythonFactory::registerType<HyperelasticForcefield_ACEgen<Element>>([template_name](sofa::core::objectmodel::Base* o) {
        return pybind11::cast(dynamic_cast<HyperelasticForcefield_ACEgen<Element>*>(o));
    });
}

template<typename Element>
void bind_hyperelastic_forcefield_ACEgen(pybind11::module &m) {
    bind_hyperelastic_forcefield_ACEgen<Element>(m, HyperelasticForcefield_ACEgen<Element>::templateName());
}

void addHyperElasticForcefield_ACEgen(pybind11::module &m);
}
