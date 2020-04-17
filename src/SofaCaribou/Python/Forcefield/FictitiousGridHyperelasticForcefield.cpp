#include "FictitiousGridHyperelasticForcefield.h"

#include <SofaCaribou/Forcefield/FictitiousGridHyperelasticForce.h>
#include <SofaCaribou/Python/Forcefield/HyperelasticForcefield.h>

namespace SofaCaribou::forcefield::python {
void addFictitiousGridHyperelasticForcefield(py::module &m) {
    using namespace SofaCaribou::forcefield;

    using caribou::geometry::SubdividedGaussHexahedron;
    using caribou::geometry::SubdividedVolumeHexahedron;

    bind_hyperelastic_forcefield<SubdividedGaussHexahedron>(m, "SubdividedGaussHexahedron");
    bind_hyperelastic_forcefield<SubdividedVolumeHexahedron>(m, "SubdividedVolumeHexahedron");

    pybind11::class_<
        FictitiousGridHyperelasticForcefield<SubdividedGaussHexahedron>,
        HyperelasticForcefield<SubdividedGaussHexahedron>,
        sofa::core::sptr<FictitiousGridHyperelasticForcefield<SubdividedGaussHexahedron>>
    > a (m, "FictitiousGridHyperelasticForcefield<SubdividedGaussHexahedron>");

    pybind11::class_<
        FictitiousGridHyperelasticForcefield<SubdividedVolumeHexahedron>,
        HyperelasticForcefield<SubdividedVolumeHexahedron>,
        sofa::core::sptr<FictitiousGridHyperelasticForcefield<SubdividedVolumeHexahedron>>
    > b (m, "FictitiousGridHyperelasticForcefield<SubdividedVolumeHexahedron>");
}
}