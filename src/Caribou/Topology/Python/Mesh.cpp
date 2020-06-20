#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include <Caribou/Constants.h>
#include <Caribou/Topology/Mesh.h>
#include <Caribou/Python/Caribou.h>
#include <Caribou/Topology/Python/Domain.h>

#include <Caribou/Geometry/Hexahedron.h>
#include <Caribou/Geometry/Quad.h>
#include <Caribou/Geometry/Segment.h>
#include <Caribou/Geometry/Triangle.h>
#include <Caribou/Geometry/Tetrahedron.h>

namespace py = pybind11;

namespace caribou::topology::python {

template <UNSIGNED_INTEGER_TYPE Dim, typename MatrixType>
void declare_mesh(py::module & m) {
    std::string name = "Mesh"+std::to_string(Dim)+"D"+typeid(MatrixType).name();
    using M = Mesh<Dim, MatrixType>;
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
//    todo(jnbrunet2000@gmail.com): Uncomment once we allow a no_copy options for the matrices (pass by ref of a np.array)
//    declare_mesh<Dim, Eigen::Matrix<float, Eigen::Dynamic, Dim, (Dim>1?Eigen::RowMajor:Eigen::ColMajor)>>(m);
    declare_mesh<Dim, Eigen::Matrix<double, Eigen::Dynamic, Dim, (Dim>1?Eigen::RowMajor:Eigen::ColMajor)>>(m);
}

void create_mesh(py::module & m) {
    declare_mesh<1>(m);
    declare_mesh<2>(m);
    declare_mesh<3>(m);
}

} // namespace caribou::topology::python