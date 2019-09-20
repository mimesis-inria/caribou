#include <pybind11/pybind11.h>

#include <pybind11/eigen.h>
#include <pybind11/stl.h>

#include <Caribou/Topology/HashGrid.h>
#include <Caribou/Geometry/Tetrahedron.h>
#include <Caribou/Geometry/Hexahedron.h>
#include <Caribou/Geometry/RectangularHexahedron.h>
#include "HashGrid.h"

namespace py = pybind11;

//PYBIND11_MAKE_OPAQUE(std::vector<py::object, std::allocator<py::object>>);

namespace caribou::topology::python
{
using namespace geometry;

template <typename Element>
auto grid(py::module & m, const std::string & name) {
    py::class_<HashGrid<Element, py::object>> o(m, name.c_str());
    o.def(py::init<FLOATING_POINT_TYPE>(), py::arg("cell_size"));
    o.def(py::init<FLOATING_POINT_TYPE, UNSIGNED_INTEGER_TYPE>(), py::arg("cell_size"), py::arg("number_of_elements"));
    o.def("add", &HashGrid<Element, py::object>::add);
    o.def("get", &HashGrid<Element, py::object>::get, py::arg("p"));
    return o;
}

void create_hashgrid(pybind11::module & m) {
    grid<Tetrahedron<interpolation::Tetrahedron4>>(m, "Tetrahedron4HashGrid");
    grid<Tetrahedron<interpolation::Tetrahedron10>>(m, "Tetrahedron10HashGrid");
    grid<Hexahedron<interpolation::Hexahedron8>>(m, "Hexahedron8HashGrid");
    grid<RectangularHexahedron<interpolation::Hexahedron8>>(m, "RectangularHexahedron8HashGrid");

    m.def("HashGrid",
        [](const py::object & ElementClass, const FLOATING_POINT_TYPE & cell_size, UNSIGNED_INTEGER_TYPE number_of_elements = 0) -> py::object {

        if (ElementClass.is(py::detail::get_type_handle(typeid(Tetrahedron<interpolation::Tetrahedron4>), true)))
            return py::cast(HashGrid<Tetrahedron<interpolation::Tetrahedron4>, py::object>(cell_size, number_of_elements));

        else if (ElementClass.is(py::detail::get_type_handle(typeid(Tetrahedron<interpolation::Tetrahedron10>), true)))
            return py::cast(HashGrid<Tetrahedron<interpolation::Tetrahedron10>, py::object>(cell_size, number_of_elements));

        else if (ElementClass.is(py::detail::get_type_handle(typeid(Hexahedron<interpolation::Hexahedron8>), true)))
            return py::cast(HashGrid<Hexahedron<interpolation::Hexahedron8>, py::object>(cell_size, number_of_elements));

        else if (ElementClass.is(py::detail::get_type_handle(typeid(RectangularHexahedron<interpolation::Hexahedron8>), true)))
            return py::cast(HashGrid<RectangularHexahedron<interpolation::Hexahedron8>, py::object>(cell_size, number_of_elements));

        return py::none();
    }, py::arg("element_type"), py::arg("cell_size"), py::arg("number_of_elements") = 0);
}

} // namespace caribou::geometry::python