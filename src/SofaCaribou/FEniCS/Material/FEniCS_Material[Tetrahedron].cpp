#include <SofaCaribou/config.h>
#include <SofaCaribou/FEniCS/Material/FEniCS_Material[Tetrahedron].h>


DISABLE_ALL_WARNINGS_BEGIN
#include <sofa/core/ObjectFactory.h>
DISABLE_ALL_WARNINGS_END

using namespace sofa::core;

using namespace caribou::geometry;

namespace SofaCaribou::material {




template class FEniCS_Material<caribou::geometry::Tetrahedron, sofa::defaulttype::Vec3Types>;
template<>
ufcx_integral* FEniCS_Material<Tetrahedron, sofa::defaulttype::Vec3Types>::FEniCS_F() {
        if(d_material_name.getValue() == "SaintVenantKirchhoff") {
            ufcx_integral *integral = form_SaintVenantKirchhoff_Tetra_F->integrals(ufcx_integral_type::cell)[0];
            return integral;
        } else if(d_material_name.getValue() == "NeoHookean") {
            ufcx_integral *integral = form_NeoHooke_Tetra_F->integrals(ufcx_integral_type::cell)[0];
            return integral;
        } else if(d_material_name.getValue() == "MooneyRivlin") {
            ufcx_integral *integral = form_MooneyRivlin_Tetra_F->integrals(ufcx_integral_type::cell)[0];
            return integral;
        } else if(d_material_name.getValue() == "Ogden") {
            ufcx_integral *integral = form_Ogden_Tetra_F->integrals(ufcx_integral_type::cell)[0];
            return integral;
        }
    }
template<>
ufcx_integral* FEniCS_Material<Tetrahedron, sofa::defaulttype::Vec3Types>::FEniCS_F_bc() {
        if(d_material_name.getValue() == "SaintVenantKirchhoff") {
            ufcx_integral *integral = form_SaintVenantKirchhoff_Tetra_F->integrals(ufcx_integral_type::exterior_facet)[0];
            return integral;
        } else if(d_material_name.getValue() == "NeoHookean") {
            ufcx_integral *integral = form_NeoHooke_Tetra_F->integrals(ufcx_integral_type::exterior_facet)[0];
            return integral;
        } else if(d_material_name.getValue() == "MooneyRivlin") {
            ufcx_integral *integral = form_MooneyRivlin_Tetra_F->integrals(ufcx_integral_type::exterior_facet)[0];
            return integral;
        } else if(d_material_name.getValue() == "Ogden") {
            ufcx_integral *integral = form_Ogden_Tetra_F->integrals(ufcx_integral_type::exterior_facet)[0];
            return integral;
        }
    }
template<>
ufcx_integral* FEniCS_Material<Tetrahedron, sofa::defaulttype::Vec3Types>::FEniCS_J() {
        if(d_material_name.getValue() == "SaintVenantKirchhoff") {
            ufcx_integral *integral = form_SaintVenantKirchhoff_Tetra_J->integrals(ufcx_integral_type::cell)[0];
            return integral;
        } else if(d_material_name.getValue() == "NeoHookean") {
            ufcx_integral *integral = form_NeoHooke_Tetra_J->integrals(ufcx_integral_type::cell)[0];
            return integral;
        } else if(d_material_name.getValue() == "MooneyRivlin") {
            ufcx_integral *integral = form_MooneyRivlin_Tetra_J->integrals(ufcx_integral_type::cell)[0];
            return integral;
        } else if(d_material_name.getValue() == "Ogden") {
            ufcx_integral *integral = form_Ogden_Tetra_J->integrals(ufcx_integral_type::cell)[0];
            return integral;
        }
    }
template<>
ufcx_integral* FEniCS_Material<Tetrahedron, sofa::defaulttype::Vec3Types>::FEniCS_Pi() {
        if(d_material_name.getValue() == "SaintVenantKirchhoff") {
            ufcx_integral *integral = form_SaintVenantKirchhoff_Tetra_Pi->integrals(ufcx_integral_type::cell)[0];
            return integral;
        } else if(d_material_name.getValue() == "NeoHookean") {
            ufcx_integral *integral = form_NeoHooke_Tetra_Pi->integrals(ufcx_integral_type::cell)[0];
            return integral;
        } else if(d_material_name.getValue() == "MooneyRivlin") {
            ufcx_integral *integral = form_MooneyRivlin_Tetra_Pi->integrals(ufcx_integral_type::cell)[0];
            return integral;
        } else if(d_material_name.getValue() == "Ogden") {
            ufcx_integral *integral = form_Ogden_Tetra_Pi->integrals(ufcx_integral_type::cell)[0];
            return integral;
        }
    }

} // namespace SofaCaribou::material

namespace sofa::core::objectmodel {
using namespace SofaCaribou::material;

[[maybe_unused]]

static int _i_ = RegisterObject("FEniCS Saint-Venant-Kirchhoff hyperelastic material")
/* Linear */         .add<FEniCS_Material<Tetrahedron, sofa::defaulttype::Vec3Types>>();

}

