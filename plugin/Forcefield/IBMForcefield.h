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
class IBMForcefield : public ForceField<DataTypes>, public DataEngine
{
public:
    SOFA_CLASS2(SOFA_TEMPLATE(IBMForcefield, DataTypes), SOFA_TEMPLATE(ForceField, DataTypes), DataEngine);

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
    using Mat6_24   = defaulttype::Mat<6, 24, Real>;
    using Mat24_6   = defaulttype::Mat<24, 6, Real>;
    using Mat24_3  = defaulttype::Mat<24, 3, Real>;
    using Mat24_24  = defaulttype::Mat<24, 24, Real>;
    using Vec3     = defaulttype::Vec<3, Real>;
    using Vec6     = defaulttype::Vec<6, Real>;
    using Vec24    = defaulttype::Vec<24, Real>;
    using PointID  =  typename HexahedronSetTopologyContainer::PointID;
    using Hexahedron   = typename HexahedronSetTopologyContainer::Hexahedron;
    using HexahedronID = size_t;
    using Triangle = sofa::core::topology::Topology::Triangle;
    using Flag = sofa::caribou::engine::CutGridEngine::Flag;

    template <typename T>
    using Link = SingleLink<IBMForcefield<DataTypes>, T, BaseLink::FLAG_STRONGLINK>;

    struct Hexa {
        std::array<PointID, 8> nodes;
        Mat24_24 K;
    };


    // Public methods
    IBMForcefield();
    void init() override;
    void reinit() override;
    void reset() override;
    void addForce(const core::MechanicalParams* mparams, Data<VecDeriv>& d_f, const Data<VecCoord>& d_x, const Data<VecDeriv>& d_v) override;
    void addDForce(const core::MechanicalParams* mparams, Data<VecDeriv>& d_df, const Data<VecDeriv>& d_dx) override;
    SReal getPotentialEnergy(const core::MechanicalParams* /* mparams */, const Data<VecCoord>& /* d_x */) const override {return 0;}
    void addKToMatrix(sofa::defaulttype::BaseMatrix * matrix, SReal kFact, unsigned int &offset) override;
    void draw(const core::visual::VisualParams* vparams) override;
    void update() override;

    static std::string templateName(const IBMForcefield<DataTypes> * = nullptr)
    {
        return DataTypes::Name();
    }

    template<class T>
    static std::string shortName(const T* ptr = nullptr, objectmodel::BaseObjectDescription* arg = nullptr)
    {
        std::string name = Inherit1::shortName(ptr, arg);
        sofa::helper::replaceAll(name, "IBMForceField", "FF");
        return name;
    }

    // Inputs
    Data< Real > d_youngModulus;
    Data< Real > d_poissonRatio;
    Data<sofa::helper::vector<Coord>> d_initial_positions;
    Data<sofa::helper::vector<Hexahedron>> d_hexahedrons;
    Data<sofa::helper::vector<Flag>> d_hexahedrons_flags;
    Data<sofa::helper::vector<Coord>> d_triangle_positions;
    Data<sofa::helper::vector<sofa::helper::vector<Triangle>>> d_triangles;

    // Output
//    Data<sofa::helper::vector<std::array<Real, 84>>> d_integrated_monomials;
    Data<sofa::helper::vector<Mat24_24>> d_stiffnesses;

private:
    std::vector<Hexa> m_hexahedrons;

};

} // namespace forcefield

} // namespace caribou

} // namespace sofa

#endif //CARIBOU_FORCEFIELD_IBMFORCEFIELD_H
