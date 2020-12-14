#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

#include <Caribou/Topology/Grid/Grid.h>

namespace py = pybind11;

namespace pybind11::detail {
using namespace caribou::topology::internal;

} // namespace pybind11::detail

namespace caribou::topology::python {
using namespace geometry;

template <class Derived>
void declare_basegrid(py::module &m) {
    // BaseGrid
    using BaseGrid = internal::BaseGrid<Derived::Dimension, Derived>;
    using WorldCoordinates = typename BaseGrid::WorldCoordinates;
    std::string name = "BaseGrid"+std::to_string(BaseGrid::Dimension)+"D";
    py::class_<BaseGrid> o (m, name.c_str());
    o.def("number_of_cells", &BaseGrid::number_of_cells);
    o.def("number_of_nodes", &BaseGrid::number_of_nodes);
    o.def("anchor_position", &BaseGrid::anchor_position);
    o.def("N", &BaseGrid::N);
    o.def("H", &BaseGrid::H);
    o.def("size", &BaseGrid::size);
    o.def("cell_index_containing", &BaseGrid::cell_index_containing, py::arg("world_coordinates"), py::arg("test_for_boundary")=true);
    o.def("cell_coordinates_at", [](const BaseGrid & g, const WorldCoordinates & c) {return g.cell_coordinates_at(c);}, py::arg("world_coordinates"));
    o.def("cell_coordinates_at", [](const BaseGrid & g, const WorldCoordinates & c, bool t) {return g.cell_coordinates_at(c, t);}, py::arg("world_coordinates"), py::arg("test_for_boundary"));
    o.def("cells_around", &BaseGrid::cells_around, py::arg("world_coordinates"));
    o.def("node", &BaseGrid::node, py::arg("node_coordinates"));
    o.def("contains", &BaseGrid::contains, py::arg("world_coordinates"), py::arg("epsilon") = EPSILON);


    // BaseDimensionalGrid is either
    //   1D    : BaseUnidimensionalGrid
    //   2D-3D :  BaseMultidimensionalGrid
    using BaseDimensionalGrid = typename Derived::Base;
    if constexpr (BaseGrid::Dimension == 1) {
        name = "BaseUnidimensionalGrid";
    } else {
        name = "BaseMultidimensionalGrid"+std::to_string(BaseGrid::Dimension)+"D";
    }

    py::class_<BaseDimensionalGrid, BaseGrid> o1 (m, name.c_str());
    o1.def("node_index_at", &BaseDimensionalGrid::node_index_at, py::arg("node_coordinates"));
    o1.def("node_coordinates_at", &BaseDimensionalGrid::node_coordinates_at, py::arg("node_index"));
    o1.def("node", &BaseDimensionalGrid::node, py::arg("node_index"));
    o1.def("cell_index_at", &BaseDimensionalGrid::cell_index_at, py::arg("cell_coordinates"));
    o1.def("cell_coordinates_at", &BaseDimensionalGrid::cell_coordinates_at, py::arg("cell_index"));
}

template <size_t Dim>
void declare_grid(py::module &m) {
    using Grid = Grid<Dim>;
    using Base = typename Grid::Base;
    using WorldCoordinates = typename Grid::WorldCoordinates;
    using GridCoordinates = typename Grid::GridCoordinates;
    using Subdivisions = typename Grid::Subdivisions;
    using Dimensions = typename Grid::Dimensions;
    using CellIndex = typename Grid::CellIndex;
    declare_basegrid<Grid>(m);
    std::string name = "Grid"+std::to_string(Grid::Dimension)+"D";
    py::class_<Grid, Base> o (m, name.c_str());
    o.def(py::init<WorldCoordinates, Subdivisions, Dimensions>(), py::arg("anchor_position"), py::arg("n"), py::arg("size"));
    o.def("node_indices_of", [](const Grid & g, const GridCoordinates & c) {return g.node_indices_of(c);}, py::arg("cell_coordinates"));
    o.def("node_indices_of", [](const Grid & g, const CellIndex & i) {return g.node_indices_of(i);}, py::arg("cell_index"));
    o.def("cell_at", [](const Grid & g, const GridCoordinates & c) {return g.cell_at(c);}, py::arg("cell_coordinates"));
    o.def("cell_at", [](const Grid & g, const CellIndex & i) {return g.cell_at(i);}, py::arg("cell_index"));

    if constexpr (Dim > 1) {
        o.def("number_of_edges", &Grid::number_of_edges);
        o.def("edge", &Grid::edge, py::arg("edge_index"));
    }

    if constexpr (Dim == 3) {
        o.def("number_of_faces", &Grid::number_of_faces);
        o.def("face", &Grid::face, py::arg("face_index"));
    }
}

void create_grid(pybind11::module & m) {
//    declare_grid<1>(m);
    declare_grid<2>(m);
    declare_grid<3>(m);
}

} // namespace caribou::topology::python