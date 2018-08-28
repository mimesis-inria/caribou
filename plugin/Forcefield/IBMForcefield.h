#ifndef CARIBOU_FORCEFIELD_IBMFORCEFIELD_H
#define CARIBOU_FORCEFIELD_IBMFORCEFIELD_H

#include <sofa/core/behavior/ForceField.h>
#include <SofaBaseTopology/HexahedronSetTopologyContainer.h>
#include "../Engine/CutGridEngine.h"

namespace sofa
{

using namespace core::objectmodel;
using namespace core::behavior;
using namespace component::topology;

namespace caribou
{

namespace forcefield
{

template<class DataTypes>
class IBMForcefield : public ForceField<DataTypes>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(IBMForcefield, DataTypes), SOFA_TEMPLATE(ForceField, DataTypes));

    // Type definitions
    using Inherit  = sofa::core::behavior::ForceField<DataTypes>;
    using VecCoord = typename DataTypes::VecCoord;
    using VecDeriv = typename DataTypes::VecDeriv;
    using Coord    = typename DataTypes::Coord;
    using Deriv    = typename DataTypes::Deriv;
    using Real     = typename Coord::value_type;
    using Mat33    = defaulttype::Mat<3, 3, Real>;
    using Mat63    = defaulttype::Mat<6, 3, Real>;
    using Mat66    = defaulttype::Mat<6, 6, Real>;
    using Mat612   = defaulttype::Mat<6, 12, Real>;
    using Mat1212  = defaulttype::Mat<12, 12, Real>;
    using Vec3     = defaulttype::Vec<3, Real>;
    using Vec6     = defaulttype::Vec<6, Real>;
    using Vec12    = defaulttype::Vec<12, Real>;
    using PointID  =  typename HexahedronSetTopologyContainer::PointID;
    using Hexahedron   = typename HexahedronSetTopologyContainer::Hexahedron;
    using HexahedronID = size_t;
    using Flag = sofa::caribou::engine::CutGridEngine::Flag;

    template <typename T>
    using Link = SingleLink<IBMForcefield<DataTypes>, T, BaseLink::FLAG_STRONGLINK>;


    // Public methods
    IBMForcefield();
    void init() override;
    void reset() override;
    void addForce(const core::MechanicalParams* mparams, Data<VecDeriv>& d_f, const Data<VecCoord>& d_x, const Data<VecDeriv>& d_v) override;
    void addDForce(const core::MechanicalParams* mparams, Data<VecDeriv>& d_df, const Data<VecDeriv>& d_dx) override;
    SReal getPotentialEnergy(const core::MechanicalParams* /* mparams */, const Data<VecCoord>& /* d_x */) const override {return 0;}
    void addKToMatrix(sofa::defaulttype::BaseMatrix * matrix, SReal kFact, unsigned int &offset) override;
    void draw(const core::visual::VisualParams* vparams) override;

    // Inputs
    Data< Real > d_youngModulus;
    Data< Real > d_poissonRatio;
    Data<sofa::helper::vector<Flag>> d_points_flags;
    Data<sofa::helper::vector<Flag>> d_hexahedrons_flags;
    Link<HexahedronSetTopologyContainer> d_container;

};

} // namespace forcefield

} // namespace caribou

} // namespace sofa

#endif //CARIBOU_FORCEFIELD_IBMFORCEFIELD_H
