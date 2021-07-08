#pragma once

#include <SofaCaribou/config.h>
#include <Caribou/Geometry/Element.h>

#include <SofaCaribou/Forcefield/CaribouForcefield.h>

DISABLE_ALL_WARNINGS_BEGIN
#include <sofa/core/behavior/MechanicalState.h>
DISABLE_ALL_WARNINGS_END

namespace SofaCaribou::forcefield {

/**
 * Defines a traction (tractive force) field.
 *
 * A traction is a force applied to a surface region (planar discretization such as triangle or quad elements). It does
 * not take into account the surface normal, but instead follows a force direction explicitly defined in the input parameters
 * of the component.
 *
 * This component allows to apply the total tractive force from a set of smaller load increments following a linear slope until
 * the total load is reach, or apply all the load at once.
 *
 * @tparam Element The element type of this force field. Will be used to determine the position
 *                 and geometric characteristics of each integration points. It must inherits from
 *                 caribou::geometry::Element.
 */
template <typename Element>
class TractionForcefield : public CaribouForcefield<Element>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(TractionForcefield, Element), SOFA_TEMPLATE(CaribouForcefield, Element));

    // Type definitions
    using Inherit  = CaribouForcefield<Element>;
    using DataTypes = typename Inherit::DataTypes;
    using VecCoord  = typename DataTypes::VecCoord;
    using VecDeriv  = typename DataTypes::VecDeriv;
    using Coord     = typename DataTypes::Coord;
    using Deriv     = typename DataTypes::Deriv;
    using Real      = typename DataTypes::Real;

    template <typename T>
    using Data = sofa::core::objectmodel::Data<T>;

    template<int nRows, int nColumns>
    using Matrix = typename caribou::geometry::Element<Element>::template Matrix<nRows, nColumns>;

    template<int nRows, int Options=0>
    using Vector = Eigen::Matrix<Real, nRows, 1, Options>;

    static constexpr INTEGER_TYPE Dimension = caribou::geometry::traits<Element>::Dimension;
    static constexpr INTEGER_TYPE NumberOfNodesPerElement = caribou::geometry::traits<Element>::NumberOfNodesAtCompileTime;
    static constexpr INTEGER_TYPE NumberOfGaussNodesPerElement = caribou::geometry::traits<Element>::NumberOfGaussNodesAtCompileTime;

    // Data structures
    struct GaussNode {
        Real weight;
        Real jacobian_determinant;
        Vector<NumberOfNodesPerElement> N;
    };

    // The container of Gauss points (for each elements) is an array if the number of integration
    // points per element is known at compile time, or a dynamic vector otherwise.
    using GaussContainer = typename std::conditional<
            NumberOfGaussNodesPerElement != caribou::Dynamic,
            std::array<GaussNode, static_cast<std::size_t>(NumberOfGaussNodesPerElement)>,
            std::vector<GaussNode>
    >::type;

    // Public methods

    CARIBOU_API
    TractionForcefield();

    CARIBOU_API
    void init() override;

    CARIBOU_API
    void reset() override;

    CARIBOU_API
    void addForce(const sofa::core::MechanicalParams* mparams, Data<VecDeriv>& d_f, const Data<VecCoord>& d_x, const Data<VecDeriv>& d_v) override;

    CARIBOU_API
    void addDForce(const sofa::core::MechanicalParams* /*mparams*/, Data<VecDeriv>& /*d_df*/, const Data<VecDeriv>& /*d_dx*/) override {}

    CARIBOU_API
    void addKToMatrix(sofa::defaulttype::BaseMatrix * /*matrix*/, SReal /*kFact*/, unsigned int & /*offset*/) override {}

    CARIBOU_API
    void handleEvent(sofa::core::objectmodel::Event* event) override;

    SReal getPotentialEnergy(const sofa::core::MechanicalParams* /*mparams*/, const Data<VecDeriv>&  /* x */) const override
    {
        msg_error() << "Get potentialEnergy not implemented";
        return 0.0;
    }

    /** Increment the traction load by an increment of traction_increment (vector of tractive force per unit area). */
    void increment_load(Deriv traction_increment_per_unit_area) ;

    /** Get the set of Gauss integration nodes of an element */
    auto gauss_nodes_of(std::size_t element_id) const -> const auto & {
        return p_elements_quadrature_nodes[element_id];
    }

private:

    // These private methods are implemented but can be overridden

    /** Compute and store the shape functions for every integration points */
    virtual void initialize_elements();

    /** Get the set of Gauss integration nodes of the given element */
    virtual auto get_gauss_nodes(const std::size_t & element_id, const Element & element) const -> GaussContainer;

private:
    // Inputs
    Data<Deriv> d_traction; ///< Tractive force per unit area
    Data<Real> d_slope; ///< Slope of force increment, the resulting traction will be p^t = p^{t-1} + p*slope. If slope = 0, the force will be constant
    Data<unsigned int> d_number_of_steps_before_increment; ///< Number of steps to wait before a load increment. This can be used to simulate a Newton-Raphson solver.

    // Outputs
    Data<VecDeriv> d_nodal_forces; /// Tractive force for each nodes of the surface on which we are applying the traction
    Data<Real> d_total_load;  ///< Current total load applied (tractive force vector times the total surface area)

    // Private variables
    std::vector<GaussContainer> p_elements_quadrature_nodes;
    bool m_traction_is_constant = false;
    Deriv m_current_traction;
    unsigned int m_number_of_steps_since_last_increment = 0;
};

} // namespace SofaCaribou::forcefield
