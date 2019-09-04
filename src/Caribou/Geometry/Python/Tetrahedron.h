#ifndef CARIBOU_GEOMETRY_PYTHON_TETRAHEDRON_H
#define CARIBOU_GEOMETRY_PYTHON_TETRAHEDRON_H

namespace pybind11{
class module;
}

namespace caribou::geometry::python
{

void create_tetrahedrons(pybind11::module & m);

} // namespace caribou::geometry::python
#endif //CARIBOU_GEOMETRY_PYTHON_TETRAHEDRON_H
