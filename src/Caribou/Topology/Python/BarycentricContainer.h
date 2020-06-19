#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

#include <Caribou/Python/Caribou.h>
#include <Caribou/Topology/BarycentricContainer.h>

namespace py = pybind11;

namespace caribou::topology::python {

template <typename Domain, typename Real>
auto interpolate_field(const BarycentricContainer<Domain> & container,
                       const Eigen::Matrix<Real, Eigen::Dynamic, Domain::Dimension> & embedded_positions,
                       const Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> & container_field_values) {
    if (container_field_values.cols() == 1) {
        // Scalar field
        Eigen::Matrix<Real, Eigen::Dynamic, 1> embedded_field_values;
        embedded_field_values.resize(container_field_values.rows(), 1);
        container.interpolate_field(embedded_positions, container_field_values, embedded_field_values);
        return py::cast(embedded_field_values);
    } else {
        Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> embedded_field_values;
        embedded_field_values.resize(container_field_values.rows(), container_field_values.cols());
        container.interpolate_field(embedded_positions, container_field_values, embedded_field_values);
        return py::cast(embedded_field_values);
    }
}

template <typename Real, typename Domain>
void declare_barycentric_container(py::class_<BarycentricContainer<Domain>> & c) {

    c.def("interpolate_field", [](const BarycentricContainer<Domain> & container,
                                  const Eigen::Matrix<Real, Eigen::Dynamic, Domain::Dimension> & embedded_positions,
                                  const Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> & container_field_values) {
        if (container_field_values.cols() == 1) {
            // Scalar field
            Eigen::Matrix<Real, Eigen::Dynamic, 1> embedded_field_values;
            embedded_field_values.resize(embedded_positions.rows(), 1);
            container.interpolate_field(embedded_positions, container_field_values, embedded_field_values);
            return py::cast(embedded_field_values);
        } else {
            Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> embedded_field_values;
            embedded_field_values.resize(embedded_positions.rows(), container_field_values.cols());
            container.interpolate_field(embedded_positions, container_field_values, embedded_field_values);
            return py::cast(embedded_field_values);
        }

    }, py::arg("embedded_positions"), py::arg("container_field_values"));
}

template <typename Domain, typename Holder>
void declare_barycentric_container(py::class_<Domain, Holder> & m) {
    std::string name = std::string("BarycentricContainer")+typeid(Domain).name();
    py::class_<BarycentricContainer<Domain>> c(m, name.c_str());

    declare_barycentric_container<float, Domain>(c);
    declare_barycentric_container<double, Domain>(c);

    m.def("BarycentricContainer", [](const Domain & domain) {
       return BarycentricContainer(&domain);
    });
}

} // namespace caribou::topology::python