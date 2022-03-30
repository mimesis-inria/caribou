#include <pybind11/pybind11.h>

#include <pybind11/eigen.h>
#include <pybind11/stl.h>

#include <Caribou/Topology/HashGrid.h>

#include <Caribou/Geometry/Tetrahedron.h>
#include <Caribou/Geometry/Tetrahedron10.h>
#include <Caribou/Geometry/Hexahedron.h>
#include <Caribou/Geometry/Hexahedron20.h>
#include <Caribou/Geometry/RectangularHexahedron.h>
#include <Caribou/Geometry/RectangularHexahedron20.h>

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
    grid<Tetrahedron>(m, "Tetrahedron4HashGrid");
    grid<Tetrahedron10>(m, "Tetrahedron10HashGrid");
    grid<Hexahedron>(m, "Hexahedron20HashGrid");
    grid<Hexahedron20>(m, "Hexahedron8HashGrid");
    grid<RectangularHexahedron>(m, "RectangularHexahedron8HashGrid");
    grid<RectangularHexahedron20>(m, "RectangularHexahedron20HashGrid");

    m.def("HashGrid",
        [](const py::object & ElementClass, const FLOATING_POINT_TYPE & cell_size, UNSIGNED_INTEGER_TYPE number_of_elements = 0) -> py::object {

        if (ElementClass.is(py::detail::get_type_handle(typeid(Tetrahedron), true)))
            return py::cast(HashGrid<Tetrahedron>(cell_size, number_of_elements));

        else if (ElementClass.is(py::detail::get_type_handle(typeid(Tetrahedron10), true)))
            return py::cast(HashGrid<Tetrahedron10>(cell_size, number_of_elements));

        else if (ElementClass.is(py::detail::get_type_handle(typeid(Hexahedron), true)))
            return py::cast(HashGrid<Hexahedron>(cell_size, number_of_elements));

        else if (ElementClass.is(py::detail::get_type_handle(typeid(Hexahedron20), true)))
            return py::cast(HashGrid<Hexahedron20>(cell_size, number_of_elements));

        else if (ElementClass.is(py::detail::get_type_handle(typeid(RectangularHexahedron), true)))
            return py::cast(HashGrid<RectangularHexahedron>(cell_size, number_of_elements));

        else if (ElementClass.is(py::detail::get_type_handle(typeid(RectangularHexahedron20), true)))
            return py::cast(HashGrid<RectangularHexahedron20>(cell_size, number_of_elements));

        return py::none();
    }, py::arg("element_type"), py::arg("cell_size"), py::arg("number_of_elements") = 0);
}

} // namespace caribou::geometry::bindings