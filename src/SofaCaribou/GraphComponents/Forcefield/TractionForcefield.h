#ifndef SOFACARIBOU_GRAPHCOMPONENTS_FORCEFIELD_TRACTIONFORCEFIELD_H
#define SOFACARIBOU_GRAPHCOMPONENTS_FORCEFIELD_TRACTIONFORCEFIELD_H

#include <sofa/defaulttype/MatSym.h>
#include <sofa/core/behavior/ForceField.h>
#include <SofaBaseTopology/TriangleSetTopologyContainer.h>
#include <sofa/core/behavior/MechanicalState.h>

namespace SofaCaribou {
namespace GraphComponents {
namespace forcefield {

using namespace sofa::core::objectmodel;
using namespace sofa::core::behavior;
using namespace sofa::component::topology;

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
 * @tparam DataTypes The datatype of the coordinates/derivatives vectors (3D float vector, 3D double vector, 2D float
 * vector or 2D double vector).
 */
template<class DataTypes>
class TractionForcefield : public sofa::core::behavior::ForceField<DataTypes>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(TractionForcefield, DataTypes), SOFA_TEMPLATE(sofa::core::behavior::ForceField, DataTypes));

    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::VecDeriv VecDeriv;
    typedef typename DataTypes::Coord    Coord   ;
    typedef typename DataTypes::Deriv    Deriv   ;
    typedef typename Coord::value_type   Real    ;
    typedef sofa::defaulttype::Mat<3,3,Real> Mat33;
    typedef sofa::defaulttype::MatSym<3,Real> MatSym3;

    using Triange = TriangleSetTopologyContainer::Triangle;

    using TriangleTopologyLink =
        SingleLink<TractionForcefield, TriangleSetTopologyContainer, BaseLink::FLAG_STRONGLINK>;

    using MechanicalStateLink =
        SingleLink<TractionForcefield, MechanicalState<DataTypes>, BaseLink::FLAG_STRONGLINK>;


    TractionForcefield();

    void init() override;
    void reset() override;
    void addForce(const sofa::core::MechanicalParams* mparams, Data<VecDeriv>& d_f, const Data<VecCoord>& d_x, const Data<VecDeriv>& d_v) override;
    void addDForce(const sofa::core::MechanicalParams* /*mparams*/, Data<VecDeriv>& /*d_df*/, const Data<VecDeriv>& /*d_dx*/) override {}
    void draw(const sofa::core::visual::VisualParams* vparams) override;
    void handleEvent(sofa::core::objectmodel::Event* event) override;
    SReal getPotentialEnergy(const sofa::core::MechanicalParams* /*mparams*/, const Data<VecDeriv>&  /* x */) const override
    {
        msg_error() << "Get potentialEnergy not implemented";
        return 0.0;
    }

    /** Increment the traction load by an increment of traction_increment (vector of tractive force per unit area). */
    void increment_load(Deriv traction_increment_per_unit_area) ;

    // Inputs
    Data<Deriv> d_traction; ///< Tractive force per unit area
    Data<sofa::helper::vector<Triange> > d_triangles; ///< List of triangles (ex: [t1p1 t1p2 t1p3 t2p1 t2p2 t2p3 ...])
    Data<Real> d_slope; ///< Slope of force increment, the resulting traction will be p^t = p^{t-1} + p*slope. If slope = 0, the force will be constant
    TriangleTopologyLink d_triangleContainer; ///< Triangle set topology container that contains the triangle indices
    MechanicalStateLink d_mechanicalState; ///< Mechanical state that contains the triangle positions
    Data<unsigned int> d_number_of_steps_before_increment; ///< Number of steps to wait before a load increment. This can be used to simulate a Newton-Raphson solver.

    // Outputs
    Data<VecDeriv> d_nodal_forces; /// Tractive force for each nodes of the surface on which we are applying the traction
    Data<Real> d_total_load;  ///< Current total load applied (tractive force vector times the total surface area)

private:

    bool m_traction_is_constant = false;
    Deriv m_current_traction;
    unsigned int m_number_of_steps_since_last_increment = 0;
};

} // namespace forcefield
} // namespace GraphComponents
} // namespace SofaCaribou


#endif //SOFACARIBOU_GRAPHCOMPONENTS_FORCEFIELD_TRACTIONFORCEFIELD_H
