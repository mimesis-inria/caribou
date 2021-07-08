#pragma once

#include <Caribou/Geometry/Element.h>
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>


namespace caribou::geometry::bindings {

template <typename Derived>
void declare_element(pybind11::module & m, const std::string & name) {
    std::string element_name = "Element"+name;
    using E = Element<Derived>;

    std::string gauss_node_name = "ElementGaussNode"+name;
    pybind11::class_<typename E::GaussNode> g(m, gauss_node_name.c_str());
    g.def_property_readonly("position", [](const typename E::GaussNode & g) {return g.position;});
    g.def_property_readonly("weight", [](const typename E::GaussNode & g) {return g.weight;});
    g.def("__str__", [](const typename E::GaussNode & g) {
        std::stringstream ss;
        Eigen::IOFormat f(std::numeric_limits<double>::digits10 + 2, 0, ", ", ", ", "[", "]", "", "");
        if (E::Dimension == 1) {
            f.matPrefix = "";
            f.matSuffix = "";
        }
        ss << "Gauss node (position = " << g.position.format(f) << ", weight = " << g.weight << ")";
        return ss.str();
    });

    pybind11::class_<E> c(m, element_name.c_str());
    c.def("number_of_nodes", &E::number_of_nodes);
    c.def("number_of_gauss_nodes", &E::number_of_gauss_nodes);
    c.def("node", &E::node);
    c.def("nodes", &E::nodes);
    c.def("gauss_node", &E::gauss_node);
    c.def("gauss_nodes", &E::gauss_nodes);
    c.def("number_of_boundary_elements", &E::number_of_boundary_elements);
    c.def("L", &E::L);
    c.def("dL", &E::dL);
    c.def("center", &E::center);
    c.def("world_coordinates", &E::world_coordinates, pybind11::arg("local_coordinates"));
    c.def("local_coordinates", [](const E & e, const typename E::WorldCoordinates & p) {return e.local_coordinates(p);}, pybind11::arg("world_coordinates"));
    c.def("contains_local", &E::contains_local, pybind11::arg("world_coordinates"), pybind11::arg("epsilon") = 1e-15);
//    c.def("interpolate", &E::interpolate);
    c.def("jacobian", &E::jacobian);
}

}