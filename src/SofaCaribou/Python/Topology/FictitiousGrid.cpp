#include <SofaCaribou/Python/Topology/FictitiousGrid.h>
#include <SofaCaribou/Topology/FictitiousGrid.inl>

#include <pybind11/stl.h>

namespace SofaCaribou::topology::python {

using sofa::defaulttype::Vec2Types;
using sofa::defaulttype::Vec3Types;

template <typename DataTypes>
auto create_fictitious_grid(py::module & m) {
    std::string name = "FictitiousGrid" + std::string(DataTypes::Name());
    py::class_<FictitiousGrid<DataTypes>, std::shared_ptr<FictitiousGrid<DataTypes>>> c (m, name.c_str());
    c.def("number_of_cells", &FictitiousGrid<DataTypes>::number_of_cells);
    c.def("number_of_nodes", &FictitiousGrid<DataTypes>::number_of_nodes);
    c.def("number_of_subdivisions", &FictitiousGrid<DataTypes>::number_of_subdivisions);
    c.def("cell_volume_ratio_distribution", &FictitiousGrid<DataTypes>::cell_volume_ratio_distribution);
    return c;
}

void addFictitiousGrid(py::module &m) {
    create_fictitious_grid<Vec2Types>(m);
    create_fictitious_grid<Vec3Types>(m);
}

} // namespace SofaCaribou::topology::python