#ifndef CARIBOU_PRESSUREFORCEFIELD_H
#define CARIBOU_PRESSUREFORCEFIELD_H

#include <sofa/defaulttype/MatSym.h>
#include <sofa/core/behavior/ForceField.h>
#include <SofaBaseTopology/TriangleSetTopologyContainer.h>
#include <sofa/core/behavior/MechanicalState.h>

namespace sofa {

namespace caribou {

namespace forcefield {

using namespace sofa::core::objectmodel;
using namespace sofa::core::behavior;
using namespace sofa::component::topology;

template<class DataTypes>
class PressureForcefield : public sofa::core::behavior::ForceField<DataTypes>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(PressureForcefield, DataTypes), SOFA_TEMPLATE(sofa::core::behavior::ForceField, DataTypes));

    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::VecDeriv VecDeriv;
    typedef typename DataTypes::Coord    Coord   ;
    typedef typename DataTypes::Deriv    Deriv   ;
    typedef typename Coord::value_type   Real    ;
    typedef sofa::defaulttype::Mat<3,3,Real> Mat33;
    typedef sofa::defaulttype::MatSym<3,Real> MatSym3;

    using Triange = TriangleSetTopologyContainer::Triangle;

    using TriangleTopologyLink =
        SingleLink<PressureForcefield, TriangleSetTopologyContainer, BaseLink::FLAG_STRONGLINK>;

    using MechanicalStateLink =
        SingleLink<PressureForcefield, MechanicalState<DataTypes>, BaseLink::FLAG_STRONGLINK>;


    PressureForcefield();

    void init() override ;
    void addForce(const sofa::core::MechanicalParams* mparams, Data<VecDeriv>& d_f, const Data<VecCoord>& d_x, const Data<VecDeriv>& d_v) override;
    void addDForce(const sofa::core::MechanicalParams* /*mparams*/, Data<VecDeriv>& /*d_df*/, const Data<VecDeriv>& /*d_dx*/) override {}
    void draw(const core::visual::VisualParams* vparams) override;
    void handleEvent(sofa::core::objectmodel::Event* event) override;
    SReal getPotentialEnergy(const core::MechanicalParams* /*mparams*/, const Data<VecDeriv>&  /* x */) const override
    {
        msg_error() << "Get potentialEnergy not implemented";
        return 0.0;
    }

    // Inputs
    Data<Deriv> d_pressure; ///< Pressure force per unit area
    Data<sofa::helper::vector<Triange> > d_triangles; ///< List of triangles (ex: [t1p1 t1p2 t1p3 t2p1 t2p2 t2p3 ...])
    Data<Real> d_slope; ///< Slope of pressure increment, the resulting pressure will be p^t = p^{t-1} + p*slope. If slope = 0, the pressure will be constant
    TriangleTopologyLink d_triangleContainer; ///< Triangle set topology container that contains the triangle indices
    MechanicalStateLink d_mechanicalState; ///< Mechanical state that contains the triangle positions
    Data<unsigned int> d_number_of_steps_before_increment; ///< Number of steps to wait before adding an increment. This can be used to simulate a Newton-Raphson solver.

    // Outputs
    Data<VecDeriv> d_nodal_forces; /// Current nodal forces from the applied pressure
    Data<Real> d_current_load;  ///< Current total load applied

private:

    void addNodalPressures(Deriv pressure) ;

    bool m_pressure_is_constant = false;
    Deriv m_current_pressure;
    unsigned int m_number_of_steps_since_last_increment = 0;
};

} // namespace forcefield

} // namespace caribou

} // namespace sofa


#endif //CARIBOU_PRESSUREFORCEFIELD_H
