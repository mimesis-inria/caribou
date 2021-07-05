#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>

#include <SofaPython3/PythonFactory.h>
#include <SofaPython3/Sofa/Core/Binding_Base.h>

#include <SofaCaribou/Mass/CaribouMass.h>
#include <SofaCaribou/Mass/CaribouMass.inl>

namespace SofaCaribou::mass::python {

template<typename Element>
void bind_caribou_mass(pybind11::module &m, const std::string & template_name) {
    pybind11::module::import("Sofa");

    std::string name = "CaribouMass<" + template_name + ">";

    pybind11::class_<CaribouMass<Element>, sofa::core::objectmodel::BaseObject, sofapython3::py_shared_ptr<CaribouMass<Element>>> c(m, name.c_str());

    c.def("M", [](const CaribouMass<Element> & self) {
        // Convert the SparseSelfAdjointView to a sparse matrix since the adjoint view isn't supported in python
        return Eigen::SparseMatrix<typename CaribouMass<Element>::Real>(self.M());
    });
    c.def("M_diag", [](const CaribouMass<Element> & self) {
        // Convert the Diagonal to a sparse matrix since the diagonal view isn't supported in python
        return Eigen::SparseMatrix<typename CaribouMass<Element>::Real>(self.M_diag());
    });
    c.def("is_diagonal", &CaribouMass<Element>::isDiagonal);
    c.def("assemble_mass_matrix", [](CaribouMass<Element> & self) {
        self.assemble_mass_matrix();
    });
    c.def("assemble_mass_matrix", [](CaribouMass<Element> & self, const Eigen::Matrix<double, Eigen::Dynamic, CaribouMass<Element>::Dimension, Eigen::RowMajor> & x) {
        self.assemble_mass_matrix(x);
    }, pybind11::arg("x").noconvert(true));
    c.def("assemble_mass_matrix", [](CaribouMass<Element> & self, const Eigen::Matrix<float, Eigen::Dynamic, CaribouMass<Element>::Dimension, Eigen::RowMajor> & x) {
        self.assemble_mass_matrix(x);
    }, pybind11::arg("x").noconvert(true));

    sofapython3::PythonFactory::registerType<CaribouMass<Element>>([template_name](sofa::core::objectmodel::Base* o) {
        return pybind11::cast(dynamic_cast<CaribouMass<Element>*>(o));
    });
}

template<typename Element>
void bind_caribou_mass(pybind11::module &m) {
    bind_caribou_mass<Element>(m, CaribouMass<Element>::templateName(static_cast<const CaribouMass<Element> *>(nullptr)));
}

void addCaribouMass(pybind11::module &m);
}