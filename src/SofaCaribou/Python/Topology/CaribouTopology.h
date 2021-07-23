#pragma once
#include <pybind11/pybind11.h>
#include <SofaCaribou/Topology/CaribouTopology.inl>
#include <SofaPython3/Sofa/Core/Binding_Base.h>

namespace SofaCaribou::topology::python {

template<typename Element>
void bind_caribou_topology(pybind11::module &m) {
    using sofa::core::objectmodel::BaseObject;
    const std::string name = "CaribouTopology<" + CaribouTopology<Element>::templateName() + ">";
    pybind11::class_<CaribouTopology<Element>, BaseObject, sofapython3::py_shared_ptr<CaribouTopology<Element>>> c(m, name.c_str());

    c.def("domain", &CaribouTopology<Element>::domain);
}

void addCaribouTopology(pybind11::module &m);

}