#pragma once

#include <SofaCaribou/config.h>

DISABLE_ALL_WARNINGS_BEGIN
#include <sofa/version.h>
#include <sofa/core/objectmodel/BaseObject.h>
#include <sofa/helper/OptionsGroup.h>
#include <sofa/core/topology/BaseTopology.h>
#include <sofa/core/behavior/ForceField.h>
DISABLE_ALL_WARNINGS_END

#include <Caribou/Geometry/Tetrahedron.h>

#if (defined(SOFA_VERSION) && SOFA_VERSION < 201200)
namespace sofa {
using Index = sofa::core::topology::Topology::index_type;
}
#endif

namespace SofaCaribou::forcefield {

using namespace sofa::core;
using namespace sofa::core::objectmodel;
using namespace sofa::core::behavior;
using namespace sofa::core::topology;
using sofa::defaulttype::Vec3Types;

class TetrahedronElasticForce : public ForceField<Vec3Types> {
public:
    SOFA_CLASS(TetrahedronElasticForce, SOFA_TEMPLATE(ForceField, Vec3Types));

    // Type definitions
    using Tetrahedron = caribou::geometry::Tetrahedron<caribou::Linear>;
    using Inherit  = ForceField<Vec3Types>;
    using DataTypes = typename Inherit::DataTypes;
    using VecCoord = typename DataTypes::VecCoord;
    using VecDeriv = typename DataTypes::VecDeriv;
    using Coord    = typename DataTypes::Coord;
    using Deriv    = typename DataTypes::Deriv;
    using Real     = typename Coord::value_type;


    template<int nRows, int nColumns, int Options=Eigen::RowMajor>
    using Matrix = Eigen::Matrix<Real, nRows, nColumns, Options>;

    template<int nRows, int nColumns, int Options=Eigen::RowMajor>
    using MatrixI = Eigen::Matrix<UNSIGNED_INTEGER_TYPE, nRows, nColumns, Options>;

    template<int nRows, int nColumns>
    using Map = Eigen::Map<const Matrix<nRows, nColumns, Eigen::RowMajor>>;

    template<int nRows, int Options=0>
    using Vector = Eigen::Matrix<Real, nRows, 1, Options>;

    template<int nRows>
    using MapVector = Eigen::Map<const Vector<nRows>>;

    using Rotation = Tetrahedron::Matrix<3,3>;
    using Mat33   =  Tetrahedron::Matrix<3,3>;
    using Vec3   =   Vector<3>;
    static constexpr INTEGER_TYPE Dimension = 3;
    static constexpr INTEGER_TYPE NumberOfNodes = Tetrahedron::NumberOfNodesAtCompileTime;

    template <typename ObjectType>
    using Link = SingleLink<TetrahedronElasticForce, ObjectType, BaseLink::FLAG_STRONGLINK>;

    // Data structures

    struct GaussNode {
        Real weight = 0;
        Real jacobian_determinant = 0;
        Matrix<NumberOfNodes, 3> dN_dx = Matrix<NumberOfNodes, 3, Eigen::RowMajor>::Zero();
    };

    // Public methods
    TetrahedronElasticForce();

    void init() override;
    void reinit() override;
    void draw(const sofa::core::visual::VisualParams* vparams) override;

    void addForce(
            const MechanicalParams* mparams,
            Data<VecDeriv>& d_f,
            const Data<VecCoord>& d_x,
            const Data<VecDeriv>& d_v) override;

    void addDForce(
            const MechanicalParams* /*mparams*/,
            Data<VecDeriv>& /*d_df*/,
            const Data<VecDeriv>& /*d_dx*/) override;

    void addKToMatrix(
            sofa::defaulttype::BaseMatrix * /*matrix*/,
            SReal /*kFact*/,
            unsigned int & /*offset*/) override;

    SReal getPotentialEnergy(
            const MechanicalParams* /* mparams */,
            const Data<VecCoord>& /* d_x */) const override
    {return 0;}

    void computeBBox(const sofa::core::ExecParams* params, bool onlyVisible) override;

private:
    /** (Re)Compute the tangent stiffness matrix */
    void compute_K();

    template <typename T>
    inline
    Tetrahedron tetrahedron(std::size_t tetrahedron_id, const T & x) const
    {
        auto * topology = d_topology_container.get();
        const auto &node_indices = topology->getTetrahedron(static_cast<sofa::Index>(tetrahedron_id));

        Matrix<NumberOfNodes, 3> m;
        for (Eigen::Index j = 0; j < NumberOfNodes; ++j) {
            const auto &node_id = node_indices[static_cast<sofa::Index>(j)];
            m.row(j) = MapVector<3>(&x[node_id][0]);
        }

        return Tetrahedron(m);
    }

private:
    Data< Real > d_youngModulus;
    Data< Real > d_poissonRatio;
    Data< bool > d_corotated;
    Link<BaseMeshTopology>   d_topology_container;

private:
    std::vector<Matrix<12, 12>> p_stiffness_matrices;
    std::vector<GaussNode> p_quadrature_nodes; // Linear tetrahedrons only have 1 gauss node per element
    std::vector<Rotation> p_initial_rotation;
    std::vector<Rotation> p_current_rotation;
};

} // namespace SofaCaribou::forcefield
