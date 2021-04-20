#pragma once

#include <SofaCaribou/config.h>

DISABLE_ALL_WARNINGS_BEGIN
#include <sofa/core/behavior/ForceField.h>
#include <sofa/core/topology/BaseTopology.h>
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/helper/OptionsGroup.h>
DISABLE_ALL_WARNINGS_END

#include <Caribou/Geometry/Hexahedron.h>

#include <Eigen/Sparse>

namespace SofaCaribou::forcefield {

using namespace sofa::core;
using namespace sofa::core::objectmodel;
using namespace sofa::core::behavior;
using namespace sofa::core::topology;
using sofa::defaulttype::Vec3Types;

class HexahedronElasticForce : public ForceField<Vec3Types>
{
public:
    SOFA_CLASS(HexahedronElasticForce, SOFA_TEMPLATE(ForceField, Vec3Types));

    // Type definitions
    using DataTypes = Vec3Types;
    using Inherit  = ForceField<DataTypes>;
    using VecCoord = typename DataTypes::VecCoord;
    using VecDeriv = typename DataTypes::VecDeriv;
    using Coord    = typename DataTypes::Coord;
    using Deriv    = typename DataTypes::Deriv;
    using Real     = typename Coord::value_type;

    using Hexahedron = caribou::geometry::Hexahedron<caribou::Linear>;
    static constexpr INTEGER_TYPE NumberOfNodes = Hexahedron::NumberOfNodesAtCompileTime;


    template<int nRows, int nColumns, int Options=Eigen::RowMajor>
    using Matrix = Eigen::Matrix<Real, nRows, nColumns, Options>;

    template<int nRows, int nColumns>
    using Map = Eigen::Map<const Matrix<nRows, nColumns, Eigen::RowMajor>>;

    template<int nRows, int Options=0>
    using Vector = Eigen::Matrix<Real, nRows, 1, Options>;

    template<int nRows>
    using MapVector = Eigen::Map<const Vector<nRows, Eigen::ColMajor>>;

    using Rotation = Hexahedron::Matrix<3,3>;
    using Mat33   = Matrix<3, 3, Eigen::RowMajor>;
    using Vec3   = Vector<3>;
    using Mat2424 = Matrix<24, 24, Eigen::RowMajor>;
    using Vec24   = Vector<24>;

    template <typename ObjectType>
    using Link = SingleLink<HexahedronElasticForce, ObjectType, BaseLink::FLAG_STRONGLINK>;

    // Data structures

    struct GaussNode {
        Real weight = 0;
        Real jacobian_determinant = 0;
        Matrix<NumberOfNodes, 3> dN_dx = Matrix<NumberOfNodes, 3>::Zero();
        Mat33 F = Mat33::Identity();
    };

    /// Integration method used to integrate the stiffness matrix.
    enum class IntegrationMethod : unsigned int {
        /// Regular 8 points gauss integration
        Regular = 0,

        /// One gauss point integration at the center of the hexahedron
        OnePointGauss = 1
    };

    // Public methods

    CARIBOU_API
    HexahedronElasticForce();

    CARIBOU_API
    void init() override;

    CARIBOU_API
    void reinit() override;

    CARIBOU_API
    void addForce(
            const MechanicalParams* mparams,
            Data<VecDeriv>& d_f,
            const Data<VecCoord>& d_x,
            const Data<VecDeriv>& d_v) override;

    CARIBOU_API
    void addDForce(
            const MechanicalParams* /*mparams*/,
            Data<VecDeriv>& /*d_df*/,
            const Data<VecDeriv>& /*d_dx*/) override;

    CARIBOU_API
    void draw(const sofa::core::visual::VisualParams* vparams) override;

    SReal getPotentialEnergy(
            const MechanicalParams* /* mparams */,
            const Data<VecCoord>& /* d_x */) const override
    {return 0;}

    CARIBOU_API
    void addKToMatrix(sofa::defaulttype::BaseMatrix * /*matrix*/, SReal /*kFact*/, unsigned int & /*offset*/) override;

    CARIBOU_API
    void computeBBox(const sofa::core::ExecParams* params, bool onlyVisible) override;

    template <typename T>
    inline
    Hexahedron hexahedron(std::size_t hexa_id, const T & x) const
    {
        auto * topology = d_topology_container.get();
        const auto &node_indices = topology->getHexahedron(static_cast<sofa::Index>(hexa_id));

        Matrix<8, 3> m;
        for (sofa::Index j = 0; j < 8; ++j) {
            const auto &node_id = node_indices[j];
            m.row(j) = MapVector<3>(&x[node_id][0]);
        }

        return Hexahedron(m);
    }

    inline
    IntegrationMethod integration_method() const
    {
        const auto m = static_cast<IntegrationMethod> (d_integration_method.getValue().getSelectedId());

        if (m == IntegrationMethod::OnePointGauss)
            return IntegrationMethod::OnePointGauss;

        return IntegrationMethod::Regular;
    }

    inline
    std::string integration_method_as_string() const
    {
        return d_integration_method.getValue().getSelectedItem();
    }

    inline
    void set_integration_method(const IntegrationMethod & m) {
        auto index = static_cast<unsigned int > (m);
        sofa::helper::WriteAccessor<Data< sofa::helper::OptionsGroup >> d = d_integration_method;
        d->setSelectedItem(index);
    }

    const std::vector<GaussNode> & gauss_nodes_of(std::size_t hexahedron_id) const {
        return p_quadrature_nodes[hexahedron_id];
    }

    const Matrix<24, 24> & stiffness_matrix_of(std::size_t hexahedron_id) const {
        return p_stiffness_matrices[hexahedron_id];
    }

    /** Get the complete tangent stiffness matrix */
    CARIBOU_API
    const Eigen::SparseMatrix<Real> & K();

    /** Get the eigen values of the tangent stiffness matrix */
    CARIBOU_API
    const Vector<Eigen::Dynamic> & eigenvalues();

    /** Get the condition number of the tangent stiffness matrix */
    CARIBOU_API
    Real cond();

private:
    /** (Re)Compute the tangent stiffness matrix */
    virtual void compute_K();

protected:
    Data< Real > d_youngModulus;
    Data< Real > d_poissonRatio;
    Data< bool > d_corotated;
    Data< sofa::helper::OptionsGroup > d_integration_method;
    Link<BaseMeshTopology>   d_topology_container;

private:
    std::vector<Matrix<24, 24>> p_stiffness_matrices;
    std::vector<std::vector<GaussNode>> p_quadrature_nodes;
    std::vector<Rotation> p_initial_rotation;
    std::vector<Rotation> p_current_rotation;
    Eigen::SparseMatrix<Real> p_K;
    Vector<Eigen::Dynamic> p_eigenvalues;
    bool K_is_up_to_date = false;
    bool eigenvalues_are_up_to_date = false;

};

} // namespace SofaCaribou::forcefield
