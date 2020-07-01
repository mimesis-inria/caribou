#include <SofaCaribou/Python/Topology/FictitiousGrid.h>
#include <SofaCaribou/Topology/FictitiousGrid.inl>

#include <SofaPython3/PythonFactory.h>
#include <SofaPython3/Sofa/Core/Binding_BaseObject.h>

#include <pybind11/stl.h>

namespace SofaCaribou::topology::python {

using sofa::defaulttype::Vec2Types;
using sofa::defaulttype::Vec3Types;

template <typename DataTypes>
auto create_fictitious_grid(py::module & m) {
    std::string name = "FictitiousGrid" + std::string(DataTypes::Name());
    py::class_<FictitiousGrid<DataTypes>, sofa::core::objectmodel::BaseObject, sofa::core::sptr<FictitiousGrid<DataTypes>>> c (m, name.c_str());
    c.def("number_of_cells", &FictitiousGrid<DataTypes>::number_of_cells);
    c.def("number_of_nodes", &FictitiousGrid<DataTypes>::number_of_nodes);
    c.def("number_of_subdivisions", &FictitiousGrid<DataTypes>::number_of_subdivisions);
    c.def("cell_volume_ratio_distribution", &FictitiousGrid<DataTypes>::cell_volume_ratio_distribution, py::arg("number_of_decimals") = 0);

    sofapython3::PythonFactory::registerType<FictitiousGrid<DataTypes>>([](sofa::core::objectmodel::Base* o) {
        return py::cast(dynamic_cast<FictitiousGrid<DataTypes>*>(o));
    });

    return c;
}

void addFictitiousGrid(py::module &m) {
    create_fictitious_grid<Vec2Types>(m);
    create_fictitious_grid<Vec3Types>(m);
}

} // namespace SofaCaribou::topology::python