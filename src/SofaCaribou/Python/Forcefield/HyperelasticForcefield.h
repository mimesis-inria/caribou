#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>

#include <SofaPython3/PythonFactory.h>
#include <SofaPython3/Sofa/Core/Binding_BaseObject.h>

#include <SofaCaribou/Forcefield/HyperelasticForcefield.h>

#include <Caribou/Constants.h>
#include <Caribou/Geometry/Tetrahedron.h>
#include <Caribou/Geometry/Hexahedron.h>

namespace py = pybind11;

template class pybind11::class_<
    SofaCaribou::forcefield::HyperelasticForcefield<caribou::geometry::Tetrahedron<caribou::Linear>>,
    sofa::core::objectmodel::BaseObject,
    sofa::core::sptr<SofaCaribou::forcefield::HyperelasticForcefield<caribou::geometry::Tetrahedron<caribou::Linear>>>
>;

template class pybind11::class_<
    SofaCaribou::forcefield::HyperelasticForcefield<caribou::geometry::Hexahedron<caribou::Linear>>,
    sofa::core::objectmodel::BaseObject,
    sofa::core::sptr<SofaCaribou::forcefield::HyperelasticForcefield<caribou::geometry::Hexahedron<caribou::Linear>>>
>;

namespace SofaCaribou::forcefield::python {

template<typename Element>
void bind_hyperelastic_forcefield(py::module &m, const std::string & template_name) {
    py::module::import("Sofa");

    std::string name = "HyperelasticForcefield<" + template_name + ">";

    py::class_<HyperelasticForcefield<Element>, sofa::core::objectmodel::BaseObject, sofa::core::sptr<HyperelasticForcefield<Element>>> c(m, name.c_str());

    c.def("K", &HyperelasticForcefield<Element>::K);
    c.def("cond", &HyperelasticForcefield<Element>::cond);
    c.def("eigenvalues", &HyperelasticForcefield<Element>::eigenvalues);

    sofapython3::PythonFactory::registerType<HyperelasticForcefield<Element>>([template_name](sofa::core::objectmodel::Base* o) {
        return py::cast(dynamic_cast<HyperelasticForcefield<Element>*>(o));
    });
}

template<typename Element>
void bind_hyperelastic_forcefield(py::module &m) {
    bind_hyperelastic_forcefield<Element>(m, HyperelasticForcefield<Element>::templateName());
}

void addHyperElasticForcefield(py::module &m);
}