#ifndef CARIBOU_TOPOLOGY_PYTHON_GRID_H
#define CARIBOU_TOPOLOGY_PYTHON_GRID_H


namespace pybind11{
class module;
}

namespace caribou::topology::python
{

void create_grid(pybind11::module & m);

} // namespace caribou::topology::python


#endif //CARIBOU_TOPOLOGY_PYTHON_GRID_H
