#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>
#include <pybind11/eigen.h>

#include <SofaPython3/bindings/Sofa/src/SofaPython3/Sofa/Core/Binding_BaseObject.h>
#include <SofaPython3/src/SofaPython3/PythonFactory.h>

#include <SofaCaribou/GraphComponents/Topology/FictitiousGrid.h>
#include <SofaCaribou/GraphComponents/Forcefield/HexahedronElasticForce.h>
#include "HexahedronElasticForce.h"

namespace py = pybind11;

using namespace sofa::core::objectmodel;

PYBIND11_MODULE(SofaCaribouPython, m) {
    m.doc() = "SofaCaribou module";

    using FictitiousGrid2D = SofaCaribou::GraphComponents::topology::FictitiousGrid<sofa::defaulttype::Vec2Types>;
    py::class_<FictitiousGrid2D, BaseObject, FictitiousGrid2D::SPtr> fictitious_grid_2d (m, "FictitiousGrid2D");
    fictitious_grid_2d.def("set_implicit_test_function", &FictitiousGrid2D::set_implicit_test_function);

    using FictitiousGrid3D = SofaCaribou::GraphComponents::topology::FictitiousGrid<sofa::defaulttype::Vec3Types>;
    py::class_<FictitiousGrid3D, BaseObject, FictitiousGrid3D::SPtr> fictitious_grid (m, "FictitiousGrid");
    fictitious_grid.def("set_implicit_test_function", &FictitiousGrid3D::set_implicit_test_function);

    SofaCaribou::Python::addHexahedronElasticForce(m);
}