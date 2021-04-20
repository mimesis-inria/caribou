#include "FictitiousGridHyperelasticForcefield.h"

#include <SofaCaribou/Forcefield/FictitiousGridHyperelasticForce.h>
#include <SofaCaribou/Python/Forcefield/HyperelasticForcefield.h>

#include <SofaPython3/Sofa/Core/Binding_Base.h>

namespace SofaCaribou::forcefield::python {
void addFictitiousGridHyperelasticForcefield(pybind11::module &m) {
    using namespace SofaCaribou::forcefield;

    using caribou::geometry::SubdividedGaussHexahedron;
    using caribou::geometry::SubdividedVolumeHexahedron;

    bind_hyperelastic_forcefield<SubdividedGaussHexahedron>(m, "SubdividedGaussHexahedron");
    bind_hyperelastic_forcefield<SubdividedVolumeHexahedron>(m, "SubdividedVolumeHexahedron");

    pybind11::class_<
        FictitiousGridHyperelasticForcefield<SubdividedGaussHexahedron>,
        HyperelasticForcefield<SubdividedGaussHexahedron>,
        sofapython3::py_shared_ptr<FictitiousGridHyperelasticForcefield<SubdividedGaussHexahedron>>
    > a (m, "FictitiousGridHyperelasticForcefield<SubdividedGaussHexahedron>");

    pybind11::class_<
        FictitiousGridHyperelasticForcefield<SubdividedVolumeHexahedron>,
        HyperelasticForcefield<SubdividedVolumeHexahedron>,
        sofapython3::py_shared_ptr<FictitiousGridHyperelasticForcefield<SubdividedVolumeHexahedron>>
    > b (m, "FictitiousGridHyperelasticForcefield<SubdividedVolumeHexahedron>");

    sofapython3::PythonFactory::registerType<FictitiousGridHyperelasticForcefield<SubdividedGaussHexahedron>>([](sofa::core::objectmodel::Base* o) {
        return pybind11::cast(dynamic_cast<FictitiousGridHyperelasticForcefield<SubdividedGaussHexahedron>*>(o));
    });

    sofapython3::PythonFactory::registerType<FictitiousGridHyperelasticForcefield<SubdividedVolumeHexahedron>>([](sofa::core::objectmodel::Base* o) {
        return pybind11::cast(dynamic_cast<FictitiousGridHyperelasticForcefield<SubdividedVolumeHexahedron>*>(o));
    });
}
}