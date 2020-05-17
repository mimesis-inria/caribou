#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include <Caribou/Topology/Mesh.h>
#include <Caribou/Python/Caribou.h>

PYBIND11_MAKE_OPAQUE(std::vector<UNSIGNED_INTEGER_TYPE>)
PYBIND11_MAKE_OPAQUE(std::vector<INTEGER_TYPE>)
PYBIND11_MAKE_OPAQUE(std::vector<FLOATING_POINT_TYPE>)
PYBIND11_MAKE_OPAQUE(std::vector<caribou::topology::Mesh<caribou::_1D>::WorldCoordinates>)
PYBIND11_MAKE_OPAQUE(std::vector<caribou::topology::Mesh<caribou::_2D>::WorldCoordinates>)
PYBIND11_MAKE_OPAQUE(std::vector<caribou::topology::Mesh<caribou::_3D>::WorldCoordinates>)

namespace py = pybind11;

namespace caribou::topology::python {

template <UNSIGNED_INTEGER_TYPE Dim>
void declare_mesh(py::module &m) {
    std::string name = "Mesh"+std::to_string(Dim)+"D";
    py::class_<Mesh<Dim>> c(m, name.c_str());

    using M = Mesh<Dim>;
    py::bind_vector<std::vector<typename M::WorldCoordinates>>(m,  "VectorWorldCoordinates"+std::to_string(Dim)+"D");

    c.def("dimension", &M::dimension);
    c.def("number_of_domains", &M::number_of_domains);
    c.def("number_of_nodes", &M::number_of_nodes);

    if constexpr (Dim == 1) {
        c.def("add_node", [](M & m, const FLOATING_POINT_TYPE & x) {
            m.add_node(typename M::WorldCoordinates(x));
        });
        c.def("add_nodes", [](M & m, const py::list & positions) {
            std::vector<typename M::WorldCoordinates> p;
            p.reserve(positions.size());
            for (const auto x : positions) {
                p.emplace_back(typename M::WorldCoordinates(py::cast<FLOATING_POINT_TYPE>(x)));
            }
            m.add_nodes(p);
        });
    }

    c.def("add_node", &M::add_node);
    c.def("add_nodes", &M::add_nodes);
    c.def("domains", &M::domains);
    c.def("domain", [](const M & m, const UNSIGNED_INTEGER_TYPE & i) {return m.domain(i);});
    c.def("domain", [](const M & m, const std::string & name) {return m.domain(name);});

    c.def("position", &M::position);
    c.def("positions", [](const M & m, const std::vector<UNSIGNED_INTEGER_TYPE> & indices) {
        return m.positions(indices);
    }, py::arg("indices").noconvert());
    c.def("positions", [](const M & m, const std::vector<INTEGER_TYPE> & indices) {
        return m.positions(indices);
    }, py::arg("indices").noconvert());
    c.def("positions", [](const M & m, const std::vector<int> & indices) {
        return m.positions(indices);
    }, py::arg("indices").noconvert());
}

void create_mesh(py::module & m) {
    py::bind_vector<std::vector<UNSIGNED_INTEGER_TYPE>>(m, "VectorULong");
    py::bind_vector<std::vector<INTEGER_TYPE>>(m, "VectorLong");
    py::bind_vector<std::vector<int>>(m, "Vectorint");
    py::bind_vector<std::vector<unsigned int>>(m, "Vectoruint");
    py::bind_vector<std::vector<FLOATING_POINT_TYPE>>(m, "VectorFloat");

    declare_mesh<1>(m);
    declare_mesh<2>(m);
    declare_mesh<3>(m);

    m.def("Mesh", [](){
        return Mesh<_3D>();
    });

    m.def("Mesh", [](const caribou::python::Dimension & dim) {
        if (dim == caribou::python::Dimension::_1D) {
            return py::cast(Mesh<_1D>());
        } else if (dim == caribou::python::Dimension::_2D) {
            return py::cast(Mesh<_2D>());
        } else {
            return py::cast(Mesh<_3D>());
        }
    });
}

} // namespace caribou::topology::python