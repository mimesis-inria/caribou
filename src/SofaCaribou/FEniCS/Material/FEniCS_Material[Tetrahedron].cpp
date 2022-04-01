#include <SofaCaribou/config.h>
#include <SofaCaribou/FEniCS/Material/FEniCS_Material[Tetrahedron].h>
//#include <SofaCaribou/Material/SaintVenantKirchhoff_Tetra.h>


DISABLE_ALL_WARNINGS_BEGIN
#include <sofa/core/ObjectFactory.h>
DISABLE_ALL_WARNINGS_END

using namespace sofa::core;

using namespace caribou::geometry;

namespace SofaCaribou::material {




template class FEniCS_Material<caribou::geometry::Tetrahedron<caribou::Linear>, sofa::defaulttype::Vec3Types>;
template<>
ufcx_integral* FEniCS_Material<Tetrahedron<caribou::Linear >, sofa::defaulttype::Vec3Types>::FEniCS_F() {
        ufcx_integral *integral = form_SaintVenantKirchhoff_Tetra_F->integrals(ufcx_integral_type::cell)[0];
        return integral;
    }
template<>
ufcx_integral* FEniCS_Material<Tetrahedron<caribou::Linear >, sofa::defaulttype::Vec3Types>::FEniCS_J() {
        ufcx_integral *integral = form_SaintVenantKirchhoff_Tetra_J->integrals(ufcx_integral_type::cell)[0];
        return integral;
    }

template class FEniCS_Material<caribou::geometry::Tetrahedron<caribou::Quadratic>, sofa::defaulttype::Vec3Types>;
template <>
ufcx_integral* FEniCS_Material<Tetrahedron<caribou::Quadratic >, sofa::defaulttype::Vec3Types>::FEniCS_F()  {
        ufcx_integral *integral = form_SaintVenantKirchhoff_Tetra_Order2_F->integrals(ufcx_integral_type::cell)[0];
        return integral;
    }
template <>
ufcx_integral* FEniCS_Material<Tetrahedron<caribou::Quadratic >, sofa::defaulttype::Vec3Types>::FEniCS_J()  {
        ufcx_integral *integral = form_SaintVenantKirchhoff_Tetra_Order2_J->integrals(ufcx_integral_type::cell)[0];
        return integral;
    }
} // namespace SofaCaribou::material

namespace sofa::core::objectmodel {
using namespace SofaCaribou::material;

[[maybe_unused]]

static int _i_ = RegisterObject("FEniCS Saint-Venant-Kirchhoff hyperelastic material")
/* Linear */         .add<FEniCS_Material<Tetrahedron<caribou::Linear >, sofa::defaulttype::Vec3Types>>()
/* Qudratic */       .add<FEniCS_Material<Tetrahedron<caribou::Quadratic >, sofa::defaulttype::Vec3Types>>();


}

