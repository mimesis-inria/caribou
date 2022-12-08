#pragma once

#include <functional>
#include <array>

#include <SofaCaribou/config.h>
#include <SofaCaribou/FEniCS/Material/FEniCS_Material.h>
#include <SofaCaribou/Forcefield/CaribouForcefield.h>

#include <Caribou/config.h>
#include <Caribou/constants.h>
#include <Caribou/Geometry/Element.h>
#include <Caribou/Topology/Mesh.h>

#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <SofaCaribou/Topology/CaribouTopology.h>
#include <sofa/core/behavior/BaseMass.h>


#if (defined(SOFA_VERSION) && SOFA_VERSION < 201200)
namespace sofa { using Index = unsigned int; }
#endif

namespace SofaCaribou::forcefield {

template <typename Element>
class HyperelasticForcefield_FEniCS : public CaribouForcefield<Element> {
public:
    SOFA_CLASS(SOFA_TEMPLATE(HyperelasticForcefield_FEniCS, Element), SOFA_TEMPLATE(CaribouForcefield, Element));

    // Type definitions
    using Inherit  = CaribouForcefield<Element>;
    using DataTypes = typename Inherit::DataTypes;
    using VecCoord  = typename DataTypes::VecCoord;
    using VecDeriv  = typename DataTypes::VecDeriv;
    using Coord     = typename DataTypes::Coord;
    using Deriv     = typename DataTypes::Deriv;
    using Real      = typename DataTypes::Real;

    static constexpr INTEGER_TYPE Dimension = caribou::geometry::traits<Element>::Dimension;
    static constexpr INTEGER_TYPE NumberOfNodesPerElement = caribou::geometry::traits<Element>::NumberOfNodesAtCompileTime;
    static constexpr INTEGER_TYPE NumberOfGaussNodesPerElement = caribou::geometry::traits<Element>::NumberOfGaussNodesAtCompileTime;

    template<int nRows, int nColumns>
    using Matrix = typename Inherit::template Matrix<nRows, nColumns>;

    template<int nRows>
    using Vector = typename Inherit::template Vector<nRows>;

    template <typename ObjectType>
    using Link = sofa::core::objectmodel::SingleLink<HyperelasticForcefield_FEniCS<Element>, ObjectType, sofa::core::objectmodel::BaseLink::FLAG_STRONGLINK>;

    
    // Public methods

    
    HyperelasticForcefield_FEniCS();

    
    void init() override;

    
    void addForce(const sofa::core::MechanicalParams* mparams, sofa::core::MultiVecDerivId fId ) override;

    
    void addForce(
        const sofa::core::MechanicalParams* mparams,
        sofa::core::objectmodel::Data<VecDeriv>& d_f,
        const sofa::core::objectmodel::Data<VecCoord>& d_x,
        const sofa::core::objectmodel::Data<VecDeriv>& d_v) override;

    
    void addDForce(
        const sofa::core::MechanicalParams* /*mparams*/,
        sofa::core::objectmodel::Data<VecDeriv>& /*d_df*/,
        const sofa::core::objectmodel::Data<VecDeriv>& /*d_dx*/) override;

    
    SReal getPotentialEnergy(
        const sofa::core::MechanicalParams* /* mparams */,
        const sofa::core::objectmodel::Data<VecCoord>& /* d_x */) const override;

    
    void addKToMatrix(sofa::defaulttype::BaseMatrix * /*matrix*/, SReal /*kFact*/, unsigned int & /*offset*/) override;


    /**
     * Get the complete tangent stiffness matrix as a compressed sparse matrix.
     *
     * \note This method will not reassembled the stiffness matrix. It will return
     *       the latest assembly done (usually during the latest Newton iteration).
     *       Use the update_stiffness() method to manually trigger a reassembly of
     *       the tangent stiffness matrix.
     * */
    auto K() const -> Eigen::SparseMatrix<Real> {
        using StorageIndex = typename Eigen::SparseMatrix<Real>::StorageIndex;

        // K is symmetric, so we only stored "one side" of the matrix.
        // But to accelerate the computation, coefficients were not
        // stored only in the upper or lower triangular part, but instead
        // in whatever triangular part (upper or lower) the first node
        // index of the element was. This means that a coefficient (i,j)
        // might be in the lower triangular part, while (k,l) is in the
        // upper triangular part. But no coefficient will be both in the
        // lower AND the upper part.

        std::vector<Eigen::Triplet<Real>> triplets;
        triplets.reserve(p_K.size()*2);
        for (StorageIndex k = 0; k < p_K.outerSize(); ++k) {
            for (typename Eigen::SparseMatrix<Real>::InnerIterator it(p_K, k); it; ++it) {
                const auto i = static_cast<StorageIndex>(it.row());
                const auto j = static_cast<StorageIndex>(it.col());
                const auto v = it.value();
                if (i != j) {
                    triplets.emplace_back(i, j, v);
                    triplets.emplace_back(j, i, v);
                } else {
                    triplets.emplace_back(i, i, v);
                }
            }
        }

        Eigen::SparseMatrix<Real> K;
        K.resize(p_K.rows(), p_K.cols());
        K.setFromTriplets(triplets.begin(), triplets.end());

        return K;
    }

    /** Get the eigen values of the tangent stiffness matrix */
    
    auto eigenvalues() -> const Vector<Eigen::Dynamic> &;

    /** Get the condition number of the tangent stiffness matrix */
    
    auto cond() -> Real;

    auto Pi() -> SReal;

    /**
     *  Assemble the stiffness matrix K.
     *
     *  Since the stiffness matrix is function of the position vector x, i.e. K(x), this method will
     *  use the mechanical state vector used in the last call to addForce (usually Position or FreePosition).
     *  If another state vector should be used as the x, use instead the update_stiffness(x) method.
     *
     *  A reference to the assembled stiffness matrix K as a column major sparse matrix can be later
     *  obtained using the method K().
     *
     */
    virtual void assemble_stiffness();

    /**
     *  Assemble the stiffness matrix K.
     *
     *  Since the stiffness matrix is function of the position vector x, i.e. K(x), this method will
     *  use the data vector x passed as parameter. If the
     *
     *  A reference to the assembled stiffness matrix K as a column major sparse matrix can be later
     *  obtained using the method K().
     *
     */
    virtual void assemble_stiffness(const sofa::core::objectmodel::Data<VecCoord> & x);

    /**
     *  Assemble the stiffness matrix K.
     *
     *  Since the stiffness matrix is function of the position vector x, i.e. K(x), this method will
     *  use the position vector x passed as a Eigen matrix nx3 parameter with n the number of nodes.
     *
     *  A reference to the assembled stiffness matrix K as a column major sparse matrix can be later
     *  obtained using the method K().
     *
     */
     template <typename Derived>
    void assemble_stiffness(const Eigen::MatrixBase<Derived> & x);

    template <typename Derived>
    void assemble_stiffness(const Eigen::MatrixBase<Derived> & x, const Eigen::MatrixBase<Derived> & x0);

private:

    // These private methods are implemented but can be overridden


    // Data members
    Link<material::FEniCS_Material<Element, DataTypes>> d_material;
    sofa::core::objectmodel::Data<bool> d_enable_multithreading;

    // Private variables
    mutable SReal Psi;
    Eigen::SparseMatrix<Real> p_K;
    Eigen::Matrix<Real, Eigen::Dynamic, 1> p_eigenvalues;

    /// Identifier of the multi-vector x used in the last call to the method addForce. This will be used to recompute
    /// the stiffness matrix K using the method update_stiffness() without any parameters.
    sofa::core::ConstMultiVecCoordId p_X_id = sofa::core::ConstVecCoordId::position();
    bool K_is_up_to_date = false;
    mutable bool Psi_is_up_to_date = false;
    bool eigenvalues_are_up_to_date = false;


    static inline double computeTotalMass(const int nb_elements, const sofa::core::behavior::BaseMass *mass) {
            double result = 0.0;
            for (int i = 0; i < nb_elements; i++) {
                result += mass->getElementMass(i);
            }
            return result;
        }
};

} // namespace SofaCaribou::forcefield
