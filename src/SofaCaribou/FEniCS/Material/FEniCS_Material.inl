#pragma once

#include <SofaCaribou/FEniCS/Material/HyperelasticMaterial_FEniCS.h>
//#include <SofaCaribou/Material/SaintVenantKirchhoff_Tetra.h>
#include <iostream>
DISABLE_ALL_WARNINGS_BEGIN
#include <SofaCaribou/FEniCS/Material/FEniCS_Generated_code/SaintVenantKirchhoff_Tetra.h>
#include <SofaCaribou/FEniCS/Material/FEniCS_Generated_code/SaintVenantKirchhoff_Tetra_Order2.h>
#include <SofaCaribou/FEniCS/Material/FEniCS_Generated_code/SaintVenantKirchhoff_Hexa.h>
#include <SofaCaribou/FEniCS/Material/FEniCS_Generated_code/SaintVenantKirchhoff_Hexa_Order2.h>

#include <SofaCaribou/FEniCS/Material/FEniCS_Generated_code/NeoHooke_Tetra.h>
#include <SofaCaribou/FEniCS/Material/FEniCS_Generated_code/NeoHooke_Tetra_Order2.h>
#include <SofaCaribou/FEniCS/Material/FEniCS_Generated_code/NeoHooke_Hexa.h>
#include <SofaCaribou/FEniCS/Material/FEniCS_Generated_code/NeoHooke_Hexa_Order2.h>
DISABLE_ALL_WARNINGS_END

namespace SofaCaribou::material {

template<typename Element, class DataTypes>
class FEniCS_Material : public HyperelasticMaterial_FEniCS< Element, DataTypes> {
    static constexpr auto Dimension = DataTypes::spatial_dimensions;
    using Coord = typename DataTypes::Coord;
    using Real  = typename Coord::value_type;
public:
    SOFA_CLASS(SOFA_TEMPLATE2(FEniCS_Material, Element, DataTypes), SOFA_TEMPLATE2(HyperelasticMaterial_FEniCS, Element, DataTypes));

    FEniCS_Material()
        : d_young_modulus(initData(&d_young_modulus,
            Real(1000), "young_modulus",
            "Young's modulus of the material",
            true /*displayed_in_GUI*/, false /*read_only_in_GUI*/))
        , d_poisson_ratio(initData(&d_poisson_ratio,
            Real(0.3),  "poisson_ratio",
           "Poisson's ratio of the material",
           true /*displayed_in_GUI*/, false /*read_only_in_GUI*/))
        , d_path(initData(&d_path,
        std::string("None"),  "path",
        "Path to FEniCS generated code",
        true /*displayed_in_GUI*/, false /*read_only_in_GUI*/))

        , d_material_name(initData(&d_material_name,
        std::string("None"),  "material_name",
        "FEniCS Material Name",
        true /*displayed_in_GUI*/, false /*read_only_in_GUI*/))
    {

    }

    Eigen::Array<float, 1, 2>
    getConstants() override {
        Eigen::Array<float, 1, 2> constants; 
        constants(0, 0) = d_young_modulus.getValue();
        constants(0, 1) = d_poisson_ratio.getValue();
        return constants;
    }

    auto FEniCS_F() -> ufcx_integral* {}
    auto FEniCS_J() -> ufcx_integral* {}
    

private:
    // Data members
    sofa::core::objectmodel::Data<Real> d_young_modulus;
    sofa::core::objectmodel::Data<Real> d_poisson_ratio;
    sofa::core::objectmodel::Data<std::string> d_path;
    sofa::core::objectmodel::Data<std::string> d_material_name;
};

} // namespace SofaCaribou::material
