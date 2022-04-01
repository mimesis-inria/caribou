#include <SofaCaribou/config.h>
#include <SofaCaribou/FEniCS/Material/SaintVenantKirchhoff/SaintVenantKirchhoffMaterial_FEniCS[Hexahedron_FEniCS].h>

DISABLE_ALL_WARNINGS_BEGIN
#include <sofa/core/ObjectFactory.h>
DISABLE_ALL_WARNINGS_END

using namespace sofa::core;

using namespace caribou::geometry;

namespace SofaCaribou::material {

template class SaintVenantKirchhoffMaterial_FEniCS<caribou::geometry::Hexahedron_FEniCS, sofa::defaulttype::Vec3Types>;

template <>
ufcx_integral* SaintVenantKirchhoffMaterial_FEniCS<Hexahedron_FEniCS, sofa::defaulttype::Vec3Types>::FEniCS_F()  {
        ufcx_integral *integral =form_SaintVenantKirchhoff_Hexa_F->integrals(ufcx_integral_type::cell)[0];
        return integral;
    }

template <>
ufcx_integral* SaintVenantKirchhoffMaterial_FEniCS<Hexahedron_FEniCS, sofa::defaulttype::Vec3Types>::FEniCS_J()  {
        ufcx_integral *integral =form_SaintVenantKirchhoff_Hexa_J->integrals(ufcx_integral_type::cell)[0];
        return integral;
    }

template <>
ufcx_integral* SaintVenantKirchhoffMaterial_FEniCS<Hexahedron_FEniCS, sofa::defaulttype::Vec3Types>::FEniCS_Pi()  {
        ufcx_integral *integral =form_SaintVenantKirchhoff_Hexa_Pi->integrals(ufcx_integral_type::cell)[0];
        return integral;
    }
} // namespace SofaCaribou::material

namespace sofa::core::objectmodel {
using namespace SofaCaribou::material;

[[maybe_unused]]

static int _i_ = RegisterObject("FEniCS Saint-Venant-Kirchhoff hyperelastic material")
/* Linear */        .add<SaintVenantKirchhoffMaterial_FEniCS<Hexahedron_FEniCS, sofa::defaulttype::Vec3Types>>();
}

