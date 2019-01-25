#ifndef SOFACARIBOU_GRAPHCOMPONENTS_FORCEFIELD_HEXAHEDRONELASTICFORCE_H
#define SOFACARIBOU_GRAPHCOMPONENTS_FORCEFIELD_HEXAHEDRONELASTICFORCE_H

#include <sofa/core/behavior/ForceField.h>
#include <sofa/core/topology/BaseTopology.h>
#include <sofa/core/behavior/MechanicalState.h>

#include <Caribou/Algebra/Matrix.h>
#include <Caribou/Algebra/Vector.h>
#include <Caribou/Geometry/Hexahedron.h>

namespace SofaCaribou::GraphComponents::forcefield {

using namespace sofa::core;
using namespace sofa::core::objectmodel;
using namespace sofa::core::behavior;
using namespace sofa::core::topology;

template<class DataTypes>
class HexahedronElasticForce : public ForceField<DataTypes>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(HexahedronElasticForce, DataTypes), SOFA_TEMPLATE(ForceField, DataTypes));

    // Type definitions
    using Inherit  = ForceField<DataTypes>;
    using VecCoord = typename DataTypes::VecCoord;
    using VecDeriv = typename DataTypes::VecDeriv;
    using Coord    = typename DataTypes::Coord;
    using Deriv    = typename DataTypes::Deriv;
    using Real     = typename Coord::value_type;

    using Mat33   = typename caribou::algebra::Matrix<3, 3, Real>;
    using Vec3   = typename caribou::algebra::Vector<3, Real>;
    using Mat2424 = typename caribou::algebra::Matrix<24, 24, Real>;
    using Vec24   = typename caribou::algebra::Vector<24, Real>;

    using Hexahedron = caribou::geometry::Hexahedron<caribou::geometry::interpolation::Hexahedron8>;

    template <typename ObjectType>
    using Link = SingleLink<HexahedronElasticForce<DataTypes>, ObjectType, BaseLink::FLAG_STRONGLINK>;

    // Public methods

    HexahedronElasticForce();

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
            const Data<VecDeriv>& /*d_dx*/) override;

    void draw(const sofa::core::visual::VisualParams* vparams) override;

    SReal getPotentialEnergy(
            const MechanicalParams* /* mparams */,
            const Data<VecCoord>& /* d_x */) const override
    {return 0;}

    void addKToMatrix(sofa::defaulttype::BaseMatrix * /*matrix*/, SReal /*kFact*/, unsigned int & /*offset*/) override {}

    void computeBBox(const sofa::core::ExecParams* params, bool onlyVisible) override;

    struct GaussNode {
        Real weight;
        Real jacobian_determinant;
        caribou::algebra::Matrix<Hexahedron::gauss_nodes.size(), 3> dN_dx;
    };

protected:
    Data< Real > d_youngModulus;
    Data< Real > d_poissonRatio;
    Data< bool > d_linear_strain;
    Link<TopologyContainer> d_topology_container;

private:
    std::vector<caribou::algebra::Matrix<24, 24, Real>> p_stiffness_matrices;
    std::vector<std::array<GaussNode,8>> p_quatrature_nodes;

};

} // namespace SofaCaribou::GraphComponents::forcefield

#endif //SOFACARIBOU_GRAPHCOMPONENTS_FORCEFIELD_HEXAHEDRONELASTICFORCE_H
