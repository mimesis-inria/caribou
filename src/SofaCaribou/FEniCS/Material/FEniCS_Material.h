#pragma once

#include <iostream>

#include <SofaCaribou/config.h>

#include <Caribou/constants.h>
#include <Caribou/Geometry/Element.h>
#include <SofaCaribou/Topology/CaribouTopology.h>
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


#include <sofa/core/objectmodel/BaseObject.h>
#include <sofa/core/behavior/MechanicalState.h>

namespace SofaCaribou::material {

template<typename Element, class DataTypes>
class FEniCS_Material : public sofa::core::objectmodel::BaseObject {
    static constexpr auto Dimension = DataTypes::spatial_dimensions;
    using Coord = typename DataTypes::Coord;
    using Real  = typename Coord::value_type;
public:
    SOFA_CLASS(SOFA_TEMPLATE2(FEniCS_Material, Element, DataTypes), sofa::core::objectmodel::BaseObject);

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


    // Retrieve FEniCS F tensor 
    virtual auto FEniCS_F() -> ufcx_integral* {}
    // Retrieve FEniCS J tensor
    virtual auto FEniCS_J() -> ufcx_integral* {}
    // Retrieve FEniCS Pi tensor
    virtual auto FEniCS_Pi() -> ufcx_integral* {}

    Eigen::Array<Real, 1, 2>
    getConstants() {
        Eigen::Array<Real, 1, 2> constants;
        constants(0, 0) = d_young_modulus.getValue();
        constants(0, 1) = d_poisson_ratio.getValue();
        return constants;
    }

    bool MaterialIsAvailable() {
        if (d_material_name.getValue() == "SaintVenantKirchhoff" || d_material_name.getValue() == "NeoHookean") {
            return true;
        } else return false;
    }

    std::string getMaterialName() {
        return d_material_name.getValue();
    }

    [[nodiscard]]  auto
    getTemplateName() const -> std::string override {
        return templateName(this);
    }

    /** Return the Element as the template name  */
    static std::string templateName(const FEniCS_Material<Element, DataTypes>* = nullptr) {
        return SofaCaribou::topology::CaribouTopology<Element>::templateName(); 
    }

    static bool canCreate(FEniCS_Material<Element, DataTypes>* o, sofa::core::objectmodel::BaseContext* context, sofa::core::objectmodel::BaseObjectDescription* arg) {
        std::string requested_data_type = arg->getAttribute( "template", "");
        std::string this_data_type = templateName(o);

        if (requested_data_type == this_data_type) {
            return true;
        }

        if (not requested_data_type.empty()) {
            arg->logError("Requested data type ('" +requested_data_type + "') is not '"+this_data_type+"'.");
            return false;
        }

        if (dynamic_cast<sofa::core::behavior::MechanicalState<DataTypes>*>(context->getMechanicalState()) != nullptr) {
            return true;
        }

        arg->logError("Cannot deduce the data type from the current context. Set the argument 'template=\""+this_data_type+"\"' to correct this.");
        return false;
    }

private:
    // Data members
    sofa::core::objectmodel::Data<Real> d_young_modulus;
    sofa::core::objectmodel::Data<Real> d_poisson_ratio;
    sofa::core::objectmodel::Data<std::string> d_path;
    sofa::core::objectmodel::Data<std::string> d_material_name;
};

} // namespace SofaCaribou::material
