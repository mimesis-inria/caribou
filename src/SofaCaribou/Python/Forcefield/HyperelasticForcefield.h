#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>

#include <SofaPython3/PythonFactory.h>
#include <SofaPython3/Sofa/Core/Binding_Base.h>

#include <SofaCaribou/Forcefield/HyperelasticForcefield.h>
#include <SofaCaribou/Forcefield/HyperelasticForcefield.inl>

namespace SofaCaribou::forcefield::python {

template<typename Element>
void bind_hyperelastic_forcefield(pybind11::module &m, const std::string & template_name) {
    pybind11::module::import("Sofa");

    std::string name = "HyperelasticForcefield<" + template_name + ">";

    pybind11::class_<HyperelasticForcefield<Element>, sofa::core::objectmodel::BaseObject, sofapython3::py_shared_ptr<HyperelasticForcefield<Element>>> c(m, name.c_str());

    c.def("K", &HyperelasticForcefield<Element>::K);
    c.def("cond", &HyperelasticForcefield<Element>::cond);
    c.def("eigenvalues", &HyperelasticForcefield<Element>::eigenvalues);
    c.def("assemble_stiffness", [](HyperelasticForcefield<Element> & self, const Eigen::Matrix<double, Eigen::Dynamic, HyperelasticForcefield<Element>::Dimension, Eigen::RowMajor> & x) {
        self.assemble_stiffness(x);
    }, pybind11::arg("x").noconvert(true));
    c.def("assemble_stiffness", [](HyperelasticForcefield<Element> & self, const Eigen::Matrix<float, Eigen::Dynamic, HyperelasticForcefield<Element>::Dimension, Eigen::RowMajor> & x) {
        self.assemble_stiffness(x);
    }, pybind11::arg("x").noconvert(true));

    sofapython3::PythonFactory::registerType<HyperelasticForcefield<Element>>([template_name](sofa::core::objectmodel::Base* o) {
        return pybind11::cast(dynamic_cast<HyperelasticForcefield<Element>*>(o));
    });
}

template<typename Element>
void bind_hyperelastic_forcefield(pybind11::module &m) {
    bind_hyperelastic_forcefield<Element>(m, HyperelasticForcefield<Element>::templateName());
}

void addHyperElasticForcefield(pybind11::module &m);
}