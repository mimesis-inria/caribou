#ifndef SOFACARIBOU_GRAPHCOMPONENTS_FORCEFIELD_TETRAHEDRONELASTICFORCE_H
#define SOFACARIBOU_GRAPHCOMPONENTS_FORCEFIELD_TETRAHEDRONELASTICFORCE_H

#include <sofa/core/objectmodel/BaseObject.h>
#include <sofa/helper/OptionsGroup.h>
#include <sofa/core/topology/BaseTopology.h>
#include <sofa/core/behavior/ForceField.h>

#include <Caribou/Geometry/Tetrahedron.h>

namespace SofaCaribou::GraphComponents::forcefield {

using namespace sofa::core;
using namespace sofa::core::objectmodel;
using namespace sofa::core::behavior;
using namespace sofa::core::topology;
using sofa::defaulttype::Vec3Types;

template<typename CanonicalTetrahedronType>
class TetrahedronElasticForce : public ForceField<Vec3Types> {
public:
    SOFA_CLASS(SOFA_TEMPLATE(TetrahedronElasticForce, CanonicalTetrahedronType), SOFA_TEMPLATE(ForceField, Vec3Types));

    // Type definitions
    using Tetrahedron = caribou::geometry::Tetrahedron<CanonicalTetrahedronType>;
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

    using Mat33   = Matrix<3, 3>;
    using Vec3   = Vector<3>;
    static constexpr INTEGER_TYPE Dimension = 3;
    static constexpr INTEGER_TYPE NumberOfNodes = Tetrahedron::NumberOfNodes;

    template <typename ObjectType>
    using Link = SingleLink<TetrahedronElasticForce<CanonicalTetrahedronType>, ObjectType, BaseLink::FLAG_STRONGLINK>;

    // Data structures

    struct GaussNode {
        Real weight = 0;
        Real jacobian_determinant = 0;
        Matrix<NumberOfNodes, 3> dN_dx = Matrix<NumberOfNodes, 3, Eigen::RowMajor>::Zero();
        Mat33 F = Mat33::Identity();
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
        const auto &node_indices = topology->getTetrahedron(tetrahedron_id);

        Matrix<NumberOfNodes, 3> m;
        for (std::size_t j = 0; j < NumberOfNodes; ++j) {
            const auto &node_id = node_indices[j];
            m.row(j) = MapVector<3>(&x[node_id][0]);
        }

        return Tetrahedron(m);
    }

private:
    Data< Real > d_youngModulus;
    Data< Real > d_poissonRatio;
    Data< bool > d_linear_strain;
    Data< bool > d_corotated;
    Link<BaseMeshTopology>   d_topology_container;

private:
    bool recompute_compute_tangent_stiffness = false;
    std::vector<Matrix<12, 12>> p_stiffness_matrices;
    std::vector<std::vector<GaussNode>> p_quadrature_nodes;
    std::vector<Mat33> p_initial_rotation;
    std::vector<Mat33> p_current_rotation;
    Eigen::SparseMatrix<Real> p_K;
    Vector<Eigen::Dynamic> p_eigenvalues;
    bool K_is_up_to_date;
    bool eigenvalues_are_up_to_date;

    Real p_mu;
    Real p_lambda;
};

} // namespace SofaCaribou::GraphComponents::forcefield

#endif //SOFACARIBOU_GRAPHCOMPONENTS_FORCEFIELD_TETRAHEDRONELASTICFORCE_H
