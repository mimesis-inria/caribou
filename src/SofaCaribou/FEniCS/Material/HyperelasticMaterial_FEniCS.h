#pragma once

#include <SofaCaribou/config.h>

#include <Caribou/constants.h>
#include <Caribou/Geometry/Element.h>
#include <SofaCaribou/Topology/CaribouTopology.h>
DISABLE_ALL_WARNINGS_BEGIN
#include <sofa/core/objectmodel/BaseObject.h>
#include <sofa/core/behavior/MechanicalState.h>

DISABLE_ALL_WARNINGS_END

#include <ufcx.h>
#include <Eigen/Eigen>
//#include <SofaCaribou/Material/SaintVenantKirchhoff_Tetra.h>


namespace SofaCaribou::material {

template<typename Element, class DataTypes>
class HyperelasticMaterial_FEniCS : public sofa::core::objectmodel::BaseObject {
    static constexpr auto Dimension = DataTypes::spatial_dimensions;
    using Coord = typename DataTypes::Coord;
    using Real  = typename Coord::value_type;
public:

    SOFA_CLASS(SOFA_TEMPLATE2(HyperelasticMaterial_FEniCS, Element, DataTypes), sofa::core::objectmodel::BaseObject);


    // Get Young Modulus and Poisson Ratio 
    virtual Eigen::Array<float, 1, 2> getConstants() {};
    // Retrieve FEniCS F tensor 
    virtual auto FEniCS_F() -> ufcx_integral* {}
    // Retrieve FEniCS J tensor
    virtual auto FEniCS_J() -> ufcx_integral* {}
    // Retrieve FEniCS Pi tensor
    virtual auto FEniCS_Pi() -> ufcx_integral* {}


    [[nodiscard]]  auto
    getTemplateName() const -> std::string override {
        return templateName(this);
    }

    /** Return the Element as the template name  */
    static std::string templateName(const HyperelasticMaterial_FEniCS<Element, DataTypes>* = nullptr) {
        return SofaCaribou::topology::CaribouTopology<Element>::templateName(); 
    }



    static bool canCreate(HyperelasticMaterial_FEniCS<Element, DataTypes>* o, sofa::core::objectmodel::BaseContext* context, sofa::core::objectmodel::BaseObjectDescription* arg) {
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

};

} // namespace SofaCaribou::material
