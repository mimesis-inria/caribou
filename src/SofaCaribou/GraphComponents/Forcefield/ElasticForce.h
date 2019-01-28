#ifndef SOFACARIBOU_GRAPHCOMPONENTS_FORCEFIELD_ELASTICFORCE_H
#define SOFACARIBOU_GRAPHCOMPONENTS_FORCEFIELD_ELASTICFORCE_H

#include <sofa/core/behavior/ForceField.h>
#include <sofa/core/topology/BaseTopology.h>
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/helper/OptionsGroup.h>

namespace SofaCaribou::GraphComponents::forcefield {

using namespace sofa::core;
using namespace sofa::core::objectmodel;
using namespace sofa::core::behavior;
using namespace sofa::core::topology;

template<class DataTypes>
class ElasticForce : public ForceField<DataTypes>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(ElasticForce, DataTypes), SOFA_TEMPLATE(ForceField, DataTypes));

    // Type definitions
    using Inherit  = ForceField<DataTypes>;
    using VecCoord = typename DataTypes::VecCoord;
    using VecDeriv = typename DataTypes::VecDeriv;
    using Coord    = typename DataTypes::Coord;
    using Deriv    = typename DataTypes::Deriv;
    using Real     = typename Coord::value_type;

    template <typename ObjectType>
    using Link = SingleLink<ElasticForce<DataTypes>, ObjectType, BaseLink::FLAG_STRONGLINK>;

    enum class ElementType {
        Edge=0, Triangle=1, Quad=2, Tetrahedron=3, Hexahedron=4, Invalid=5
    };

    // Public methods

    ElasticForce();

    void init() override;
    void reinit() override;

    void addForce(
            const MechanicalParams* mparams,
            Data<VecDeriv>& d_f,
            const Data<VecCoord>& d_x,
            const Data<VecDeriv>& d_v) override;

    void addDForce(
            const MechanicalParams* /*mparams*/,
            Data<VecDeriv>& /*d_df*/,
            const Data<VecDeriv>& /*d_dx*/) override {}

    void draw(const sofa::core::visual::VisualParams* vparams) override;

    SReal getPotentialEnergy(
            const MechanicalParams* /* mparams */,
            const Data<VecCoord>& /* d_x */) const override
    {return 0;}

    void addKToMatrix(sofa::defaulttype::BaseMatrix * /*matrix*/, SReal /*kFact*/, unsigned int & /*offset*/) override {}

    void computeBBox(const sofa::core::ExecParams* params, bool onlyVisible) override;

    inline std::string element_type_string() const {
        switch (m_element_type) {
            case ElementType::Edge:        return "Edge";
            case ElementType::Triangle:    return "Triangle";
            case ElementType::Quad:        return "Quad";
            case ElementType::Tetrahedron: return "Tetrahedron";
            case ElementType::Hexahedron:  return "Hexahedron";
            default: return "Unknown";
        }
    }

    inline ElementType element_type(size_t option_index) {
        switch (option_index) {
            case 0: return ElementType::Edge;
            case 1: return ElementType::Triangle;
            case 2: return ElementType::Quad;
            case 3: return ElementType::Tetrahedron;
            case 4: return ElementType::Hexahedron;
            default:return ElementType::Invalid;
        }
    }

protected:
    Data< Real > d_youngModulus;
    Data< Real > d_poissonRatio;
    Data<sofa::helper::OptionsGroup> d_element_type;
    Link<TopologyContainer> d_topology_container;

private:
    ElementType m_element_type;

};

} // namespace SofaCaribou::GraphComponents::forcefield

#endif //SOFACARIBOU_GRAPHCOMPONENTS_FORCEFIELD_ELASTICFORCE_H
