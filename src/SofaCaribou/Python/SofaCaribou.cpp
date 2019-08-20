#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>

#include <SofaPython3/bindings/Sofa/src/SofaPython3/Sofa/Core/Binding_BaseObject.h>

#include <SofaCaribou/GraphComponents/Topology/FictitiousGrid.h>
#include <SofaCaribou/GraphComponents/Forcefield/HexahedronElasticForce.h>
#include "HexahedronElasticForce.h"

namespace py = pybind11;

using namespace sofa::core::objectmodel;

PYBIND11_MODULE(SofaCaribouPython, m) {
    m.doc() = "SofaCaribou module";

    using FictitiousGrid3D = SofaCaribou::GraphComponents::topology::FictitiousGrid<sofa::defaulttype::Vec3Types>;
    py::class_<FictitiousGrid3D> fictitious_grid (m, "FictitiousGrid");
    fictitious_grid.def("set_implicit_test_function", &FictitiousGrid3D::set_implicit_test_function);

    SofaCaribou::Python::addHexahedronElasticForce(m);
}