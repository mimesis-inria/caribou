#ifndef CARIBOU_GEOMETRY_PYTHON_POINT_H
#define CARIBOU_GEOMETRY_PYTHON_POINT_H

namespace pybind11{
class module;
}

namespace caribou
{
namespace geometry
{
namespace python
{

void create_point(pybind11::module & m);

} // namespace python

} // namespace geometry

} // namespace caribou
#endif //CARIBOU_GEOMETRY_PYTHON_POINT_H
