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

#include <SofaCaribou/FEniCS/Material/FEniCS_Generated_code/MooneyRivlin_Tetra.h>
#include <SofaCaribou/FEniCS/Material/FEniCS_Generated_code/MooneyRivlin_Tetra_Order2.h>
#include <SofaCaribou/FEniCS/Material/FEniCS_Generated_code/MooneyRivlin_Hexa.h>
#include <SofaCaribou/FEniCS/Material/FEniCS_Generated_code/MooneyRivlin_Hexa_Order2.h>

#include <SofaCaribou/FEniCS/Material/FEniCS_Generated_code/Ogden_Tetra.h>
#include <SofaCaribou/FEniCS/Material/FEniCS_Generated_code/Ogden_Tetra_Order2.h>
#include <SofaCaribou/FEniCS/Material/FEniCS_Generated_code/Ogden_Hexa.h>
#include <SofaCaribou/FEniCS/Material/FEniCS_Generated_code/Ogden_Hexa_Order2.h>
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
        , d_C01(initData(&d_C01, 
            Real(0), "C01", 
            "Monoey-Rivlin parameter c1", 
            true /*displayed_in_GUI*/, false /*read_only_in_GUI*/))
        , d_C10(initData(&d_C10, 
            Real(0), "C10", 
            "Mooney-Rivlin parameter c2", 
            true /*displayed_in_GUI*/, false /*read_only_in_GUI*/))
        , d_k(initData(&d_k,
            Real(0), "k",
            "Mooney-Rivlin parameter k",
            true /*displayed_in_GUI*/, false /*read_only_in_GUI*/))
        , d_bulk_modulus(initData(&d_bulk_modulus,
            Real(0), "bulk_modulus",
            "Ogden parameter bulk modulus",
            true /*displayed_in_GUI*/, false /*read_only_in_GUI*/))
        , d_a(initData(&d_a,
            Real(0), "a",
            "Ogden parameter a",
            true /*displayed_in_GUI*/, false /*read_only_in_GUI*/))
        , d_b(initData(&d_b,
            Real(0), "b",
            "Ogden parameter b",
            true /*displayed_in_GUI*/, false /*read_only_in_GUI*/))
        , d_a_f(initData(&d_a_f,
            Real(0), "a_f",
            "Ogden parameter a_f",
            true /*displayed_in_GUI*/, false /*read_only_in_GUI*/))
        , d_b_f(initData(&d_b_f,
            Real(0), "b_f",
            "Ogden parameter b_f",
            true /*displayed_in_GUI*/, false /*read_only_in_GUI*/))
        , d_a_s(initData(&d_a_s,
            Real(0), "a_s",
            "Ogden parameter a_s",
            true /*displayed_in_GUI*/, false /*read_only_in_GUI*/))
        , d_b_s(initData(&d_b_s,
            Real(0), "b_s",
            "Ogden parameter b_s",
            true /*displayed_in_GUI*/, false /*read_only_in_GUI*/))
        , d_a_fs(initData(&d_a_fs,
            Real(0), "a_fs",
            "Ogden parameter a_fs",
            true /*displayed_in_GUI*/, false /*read_only_in_GUI*/))
        , d_b_fs(initData(&d_b_fs,
            Real(0), "b_fs",
            "Ogden parameter b_fs",
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
    getYoungModulusAndPoissonRatio() {
        Eigen::Array<Real, 1, 2> constants;
        constants(0, 0) = d_young_modulus.getValue();
        constants(0, 1) = d_poisson_ratio.getValue();
        return constants;
    }


    Eigen::Array<Real, 1, 3> getMooneyRivlinConstants() {
        Eigen::Array<Real, 1, 3> constants;
        constants(0, 0) = d_C01.getValue();
        constants(0, 1) = d_C10.getValue();
        constants(0, 2) = d_k.getValue();
        return constants;
    }

    Eigen::Array<Real, 1, 9> getOgdenConstants() {
        Eigen::Array<Real, 1, 9> constants;
        constants(0, 0) = d_bulk_modulus.getValue();
        constants(0, 1) = d_a.getValue();
        constants(0, 2) = d_b.getValue();
        constants(0, 3) = d_a_f.getValue();
        constants(0, 4) = d_b_f.getValue();
        constants(0, 5) = d_a_s.getValue();
        constants(0, 6) = d_b_s.getValue();
        constants(0, 7) = d_a_fs.getValue();
        constants(0, 8) = d_b_fs.getValue();
        return constants;
    }


    bool MaterialIsAvailable() {
        if (d_material_name.getValue() == "SaintVenantKirchhoff" || 
            d_material_name.getValue() == "NeoHookean" ||
            d_material_name.getValue() == "MooneyRivlin" ||
            d_material_name.getValue() == "Ogden") {
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
    sofa::core::objectmodel::Data<Real> d_C01;
    sofa::core::objectmodel::Data<Real> d_C10;
    sofa::core::objectmodel::Data<Real> d_k;
    sofa::core::objectmodel::Data<Real> d_bulk_modulus;
    sofa::core::objectmodel::Data<Real> d_a;
    sofa::core::objectmodel::Data<Real> d_b;
    sofa::core::objectmodel::Data<Real> d_a_f;
    sofa::core::objectmodel::Data<Real> d_b_f;
    sofa::core::objectmodel::Data<Real> d_a_s;
    sofa::core::objectmodel::Data<Real> d_b_s;
    sofa::core::objectmodel::Data<Real> d_a_fs;
    sofa::core::objectmodel::Data<Real> d_b_fs;
    sofa::core::objectmodel::Data<std::string> d_path;
    sofa::core::objectmodel::Data<std::string> d_material_name;
};

} // namespace SofaCaribou::material
