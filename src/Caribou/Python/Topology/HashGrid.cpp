#include <pybind11/pybind11.h>

#include <pybind11/eigen.h>
#include <pybind11/stl.h>

#include <Caribou/Topology/HashGrid.h>

#include <Caribou/Geometry/Tetrahedron.h>
#include <Caribou/Geometry/Hexahedron.h>
#include <Caribou/Geometry/RectangularHexahedron.h>

namespace py = pybind11;

//PYBIND11_MAKE_OPAQUE(std::vector<py::object, std::allocator<py::object>>);

namespace caribou::topology::bindings
{
using namespace geometry;

template <typename Element>
auto grid(py::module & m, const std::string & name) {
    py::class_<HashGrid<Element>> o(m, name.c_str());
    o.def(py::init<FLOATING_POINT_TYPE>(), py::arg("cell_size"));
    o.def(py::init<FLOATING_POINT_TYPE, UNSIGNED_INTEGER_TYPE>(), py::arg("cell_size"), py::arg("number_of_elements"));
    o.def("add", &HashGrid<Element>::add);
    o.def("get", &HashGrid<Element>::get, py::arg("p"));
    return o;
}

void create_hashgrid(pybind11::module & m) {
    grid<Tetrahedron<Linear>>(m, "Tetrahedron4HashGrid");
    grid<Tetrahedron<Quadratic>>(m, "Tetrahedron10HashGrid");
    grid<Hexahedron<Linear>>(m, "Hexahedron20HashGrid");
    grid<Hexahedron<Quadratic>>(m, "Hexahedron8HashGrid");
    grid<RectangularHexahedron<Linear>>(m, "RectangularHexahedron8HashGrid");
    grid<RectangularHexahedron<Quadratic>>(m, "RectangularHexahedron20HashGrid");

    m.def("HashGrid",
        [](const py::object & ElementClass, const FLOATING_POINT_TYPE & cell_size, UNSIGNED_INTEGER_TYPE number_of_elements = 0) -> py::object {

        if (ElementClass.is(py::detail::get_type_handle(typeid(Tetrahedron<Linear>), true)))
            return py::cast(HashGrid<Tetrahedron<Linear>>(cell_size, number_of_elements));

        else if (ElementClass.is(py::detail::get_type_handle(typeid(Tetrahedron<Quadratic>), true)))
            return py::cast(HashGrid<Tetrahedron<Quadratic>>(cell_size, number_of_elements));

        else if (ElementClass.is(py::detail::get_type_handle(typeid(Hexahedron<Linear>), true)))
            return py::cast(HashGrid<Hexahedron<Linear>>(cell_size, number_of_elements));

        else if (ElementClass.is(py::detail::get_type_handle(typeid(Hexahedron<Quadratic>), true)))
            return py::cast(HashGrid<Hexahedron<Quadratic>>(cell_size, number_of_elements));

        else if (ElementClass.is(py::detail::get_type_handle(typeid(RectangularHexahedron<Linear>), true)))
            return py::cast(HashGrid<RectangularHexahedron<Linear>>(cell_size, number_of_elements));

        else if (ElementClass.is(py::detail::get_type_handle(typeid(RectangularHexahedron<Quadratic>), true)))
            return py::cast(HashGrid<RectangularHexahedron<Quadratic>>(cell_size, number_of_elements));

        return py::none();
    }, py::arg("element_type"), py::arg("cell_size"), py::arg("number_of_elements") = 0);
}

} // namespace caribou::geometry::bindings