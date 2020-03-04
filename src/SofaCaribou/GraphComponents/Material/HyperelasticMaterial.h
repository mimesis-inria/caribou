#pragma once

#include <Caribou/config.h>
#include <Eigen/Eigen>
#include <sofa/core/objectmodel/BaseObject.h>
#include <sofa/core/behavior/MechanicalState.h>

namespace SofaCaribou::GraphComponents::material {

template<class DataTypes>
class HyperelasticMaterial : public sofa::core::objectmodel::BaseObject {
    static constexpr auto Dimension = DataTypes::spatial_dimensions;
    using Coord = typename DataTypes::Coord;
    using Real  = typename Coord::value_type;
public:

    SOFA_CLASS(SOFA_TEMPLATE(HyperelasticMaterial, DataTypes), sofa::core::objectmodel::BaseObject);

    /**
     * This is called just before the material is updated on every points (usually just before a Newton step).
     * It can be used to update some coefficients that will be used on material points (for example, compute
     * the parameters mu and lambda from the young modulus and poisson ratio given as data arguments).
     */
    virtual void
    before_update() {}

    /**
     * Get the strain energy density Psi from the Green-Lagrange strain tensor E.
     */
    virtual Real
    strain_energy_density(const Real & J, const Eigen::Matrix<Real, Dimension, Dimension>  & E) const = 0;

    /**
     * Get the second Piola-Kirchhoff stress tensor from the Green-Lagrange strain tensor E.
     *
     * With the energy density function Psi(E) -> scalar
     * the second Piola-Kirchhoff stress tensor S is defined as
     *   S = d(Psi)/dE
     * and is a symmetric second order matrix (eg. 3x3 matrix).
     */
    virtual Eigen::Matrix<Real, Dimension, Dimension>
    PK2_stress(const Real & J, const Eigen::Matrix<Real, Dimension, Dimension>  & E) const = 0;

    /**
     * Get the jacobian of the second Piola-Kirchhoff stress tensor w.r.t the Green-Lagrange strain tensor E.
     *
     * With the energy density function Psi(E) -> scalar
     * the jacobian of the second Piola-Kirchhoff stress tensor S is defined as
     *   J = d^2(Psi)/dE^2 = d(S)/dE
     * and is a symmetric fourth order matrix represented in its compressed format (eg. 6x6 matrix).
     */
    virtual Eigen::Matrix<Real, 6, 6>
    PK2_stress_jacobian(const Real & J, const Eigen::Matrix<Real, Dimension, Dimension>  & E) const = 0;


    // Sofa's scene methods

    /** Return the data type (ex. Vec3D) as the template name  */
    static std::string templateName(const HyperelasticMaterial<DataTypes>* = nullptr) {
        return DataTypes::Name();
    }

    /**
     * Check if we can create a HyperelasticMaterial<DataTypes> with the given
     * template argument (from the scene parser).
     *
     * If no template argument was specified to the component, try to find a
     * mechanical state in the current context to deduce the data type.
     */
    static bool canCreate(HyperelasticMaterial<DataTypes>* o, sofa::core::objectmodel::BaseContext* context, sofa::core::objectmodel::BaseObjectDescription* arg) {
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

} // namespace SofaCaribou::GraphComponents::material
