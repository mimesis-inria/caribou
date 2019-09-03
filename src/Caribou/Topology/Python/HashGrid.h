#ifndef CARIBOU_TOPOLOGY_PYTHON_HASHGRID_H
#define CARIBOU_TOPOLOGY_PYTHON_HASHGRID_H

namespace pybind11{
class module;
}

namespace caribou::topology::python
{

void create_hashgrid(pybind11::module & m);

} // namespace caribou::topology::python
#endif //CARIBOU_TOPOLOGY_PYTHON_HASHGRID_H
