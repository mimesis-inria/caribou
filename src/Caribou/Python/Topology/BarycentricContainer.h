#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

#include <Caribou/Topology/BarycentricContainer.h>
#include <Caribou/Python/Caribou.h>

namespace py = pybind11;

namespace caribou::topology::bindings {

template <typename Real, typename Domain>
void declare_barycentric_container(py::class_<BarycentricContainer<Domain>> & c) {

    c.def("interpolate", [](const BarycentricContainer<Domain> & container,
                            const Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> & container_field_values) {
        if (container_field_values.cols() == 1) {
            // Scalar field
            Eigen::Matrix<Real, Eigen::Dynamic, 1> embedded_field_values;
            embedded_field_values.resize(container.barycentric_points().size(), 1);
            container.interpolate(container_field_values, embedded_field_values);
            return py::cast(embedded_field_values);
        } else {
            Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> embedded_field_values;
            embedded_field_values.resize(container.barycentric_points().size(), container_field_values.cols());
            container.interpolate(container_field_values, embedded_field_values);
            return py::cast(embedded_field_values);
        }

    }, py::arg("container_field_values"));

    c.def_property_readonly("outside_nodes", [](const BarycentricContainer<Domain> & container){
        return container.outside_nodes();
    });

    c.def_property_readonly("barycentric_points", [](const BarycentricContainer<Domain> & container){
        return container.barycentric_points();
    });

    c.def("closest_elements", [](const BarycentricContainer<Domain> & container, const typename BarycentricContainer<Domain>::WorldCoordinates & p){
        return container.closest_elements(p);
    }, py::arg("world_coordinates"));
}

template <typename Domain, typename Holder>
void declare_barycentric_container(py::class_<Domain, Holder> & m) {
    std::string name = std::string("BarycentricContainer")+typeid(Domain).name();
    py::class_<BarycentricContainer<Domain>> c(m, name.c_str());

    // BarycentricPoint
    std::string bp_name = name+"_BarycentricPoint";
    py::class_<typename BarycentricContainer<Domain>::BarycentricPoint> bp(c, name.c_str());
    bp.def_property_readonly("element_index", [](const typename BarycentricContainer<Domain>::BarycentricPoint & bp) {
        return bp.element_index;
    });
    bp.def_property_readonly("local_coordinates", [](const typename BarycentricContainer<Domain>::BarycentricPoint & bp) {
        return bp.local_coordinates;
    });

//    todo(jnbrunet2000@gmail.com): Uncomment once we allow a no_copy options for the matrices (pass by ref of a np.array)
//    declare_barycentric_container<float, Domain>(c);
    declare_barycentric_container<double, Domain>(c);

    m.def("embed", [](const Domain & domain, const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> & embedded_points) {
       return domain.embed(embedded_points);
    }, py::arg("embedded_points"));
}

} // namespace caribou::topology::bindings