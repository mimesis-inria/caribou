#ifndef CARIBOU_GEOMETRY_PYTHON_HEXAHEDRON_H
#define CARIBOU_GEOMETRY_PYTHON_HEXAHEDRON_H

namespace pybind11{
class module;
}

namespace caribou::geometry::python
{

void create_hexahedrons(pybind11::module & m);

} // namespace caribou::geometry::python
#endif //CARIBOU_GEOMETRY_PYTHON_HEXAHEDRON_H
