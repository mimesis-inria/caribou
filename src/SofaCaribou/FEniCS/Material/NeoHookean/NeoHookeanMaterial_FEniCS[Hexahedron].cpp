#include <SofaCaribou/config.h>
#include <SofaCaribou/FEniCS/Material/NeoHookean/NeoHookeanMaterial_FEniCS[Hexahedron].h>

DISABLE_ALL_WARNINGS_BEGIN
#include <sofa/core/ObjectFactory.h>
DISABLE_ALL_WARNINGS_END

using namespace sofa::core;

using namespace caribou::geometry;

namespace SofaCaribou::material {

template class NeoHookeanMaterial_FEniCS<caribou::geometry::Hexahedron<caribou::Linear>, sofa::defaulttype::Vec3Types>;

template <>
ufcx_integral* NeoHookeanMaterial_FEniCS<Hexahedron<caribou::Linear >, sofa::defaulttype::Vec3Types>::FEniCS_F()  {
        ufcx_integral *integral =form_NeoHooke_Hexa_F->integrals(ufcx_integral_type::cell)[0];
        return integral;
    }
template <>
ufcx_integral* NeoHookeanMaterial_FEniCS<Hexahedron<caribou::Linear >, sofa::defaulttype::Vec3Types>::FEniCS_J()  {
        ufcx_integral *integral =form_NeoHooke_Hexa_J->integrals(ufcx_integral_type::cell)[0];
        return integral;
    }
template <>
ufcx_integral* NeoHookeanMaterial_FEniCS<Hexahedron<caribou::Linear >, sofa::defaulttype::Vec3Types>::FEniCS_Pi()  {
        ufcx_integral *integral =form_NeoHooke_Hexa_Pi->integrals(ufcx_integral_type::cell)[0];
        return integral;
    }



template class NeoHookeanMaterial_FEniCS<caribou::geometry::Hexahedron<caribou::Quadratic>, sofa::defaulttype::Vec3Types>;
template <>
ufcx_integral* NeoHookeanMaterial_FEniCS<Hexahedron<caribou::Quadratic >, sofa::defaulttype::Vec3Types>::FEniCS_F()  {
        ufcx_integral *integral = form_NeoHooke_Hexa_Order2_F->integrals(ufcx_integral_type::cell)[0];
        return integral;
    }
template <>
ufcx_integral* NeoHookeanMaterial_FEniCS<Hexahedron<caribou::Quadratic >, sofa::defaulttype::Vec3Types>::FEniCS_J()  {
        ufcx_integral *integral = form_NeoHooke_Hexa_Order2_J->integrals(ufcx_integral_type::cell)[0];
        return integral;
    }
template <>
ufcx_integral* NeoHookeanMaterial_FEniCS<Hexahedron<caribou::Quadratic >, sofa::defaulttype::Vec3Types>::FEniCS_Pi()  {
        ufcx_integral *integral = form_NeoHooke_Hexa_Order2_Pi->integrals(ufcx_integral_type::cell)[0];
        return integral;
    }


} // namespace SofaCaribou::material

namespace sofa::core::objectmodel {
using namespace SofaCaribou::material;

[[maybe_unused]]

static int _i_ = RegisterObject("FEniCS NeoHookean hyperelastic material")
/* Linear */        .add<NeoHookeanMaterial_FEniCS<Hexahedron<caribou::Linear >, sofa::defaulttype::Vec3Types>>()
/* Quadratic*/      .add<NeoHookeanMaterial_FEniCS<Hexahedron<caribou::Quadratic >, sofa::defaulttype::Vec3Types>>();

}

