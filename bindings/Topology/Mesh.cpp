#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include <Caribou/constants.h>
#include <Caribou/Topology/Mesh.h>
#include <Caribou/Bindings/Topology/Domain.h>

namespace py = pybind11;

namespace caribou::topology::bindings {

template <UNSIGNED_INTEGER_TYPE Dim, typename MatrixType>
void declare_mesh(py::module & m) {
    using M = Mesh<Dim, EigenNodesHolder<MatrixType>>;
    std::string name = typeid(M).name();
    py::class_<M> c(m, name.c_str());

    c.def("dimension", &M::dimension);
    c.def("number_of_domains", &M::number_of_domains);
    c.def("number_of_nodes", &M::number_of_nodes);

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

    declare_domains(c);

    m.def("Mesh", [](const MatrixType & nodes) {
        return py::cast(M(nodes));
    }, py::arg("nodes"));
}

template <UNSIGNED_INTEGER_TYPE Dim>
void declare_mesh(py::module &m) {
//    declare_mesh<Dim, Eigen::Matrix<float, Eigen::Dynamic, Dim, (Dim>1?Eigen::RowMajor:Eigen::ColMajor)>>(m);
    declare_mesh<Dim, Eigen::Matrix<double, Eigen::Dynamic, Dim, (Dim>1?Eigen::RowMajor:Eigen::ColMajor)>>(m);
}

void create_mesh(py::module & m) {
    declare_mesh<1>(m);
    declare_mesh<2>(m);
    declare_mesh<3>(m);
}

} // namespace caribou::topology::bindings