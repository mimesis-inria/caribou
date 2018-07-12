#ifndef CARIBOU_FORCEFIELD_FEMFORCEFIELD_H
#define CARIBOU_FORCEFIELD_FEMFORCEFIELD_H

#include <sofa/core/behavior/ForceField.h>
#include <sofa/defaulttype/Vec.h>
#include <sofa/helper/OptionsGroup.h>
#include <sofa/helper/decompose.h>
#include <SofaBaseTopology/TetrahedronSetTopologyContainer.h>

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
class FEMForcefield : public ForceField<DataTypes>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(FEMForcefield, DataTypes), SOFA_TEMPLATE(ForceField, DataTypes));

    // Enumerations
    enum class Corotational {
        NONE = 0,    ///< No corotational method used (small displacements)
        LARGE = 1,   ///< Corotational method based on a QR decomposition    -> Nesme et al 2005 "Efficient, Physically Plausible Finite Elements"
        POLAR = 2,   ///< Corotational method based on a polar decomposition -> Muller et al 2004 "Interactive Virtual Materials"
        SVD = 3,      ///< Corotational method based on a SVD decomposition   -> inspired from Irving et al 2004 "Invertible Finite Element for Robust Simulation of Large Deformation"
        DEFORMATION_GRADIENT = 4  ///< Corotational method based on rotation extraction from the deformation gradient (F = GradU + I)
    };

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
    using PointID  =  typename TetrahedronSetTopologyContainer::PointID;
    using Tetrahedron   = typename TetrahedronSetTopologyContainer::Tetrahedron;
    using TetrahedronID = size_t;


    // Public methods
    FEMForcefield();
    void init() override;
    void addForce(const core::MechanicalParams* mparams, Data<VecDeriv>& d_f, const Data<VecCoord>& d_x, const Data<VecDeriv>& d_v) override;
    void addDForce(const core::MechanicalParams* mparams, Data<VecDeriv>& d_df, const Data<VecDeriv>& d_dx) override;
    SReal getPotentialEnergy(const core::MechanicalParams* /* mparams */, const Data<VecCoord>& /* d_x */) const override {return 0;}
    void addKToMatrix(sofa::defaulttype::BaseMatrix * matrix, SReal kFact, unsigned int &offset) override;

    inline Mat612 getStrainDisplacement(const helper::fixed_array<Deriv, 4> shapeDerivatives) const
    {
        Mat612 B;

        for (unsigned int i = 0; i < 4; ++i) {
            const Deriv & GradS = shapeDerivatives[i];

            B[0][i*3+0] = GradS[0];  B[0][i*3+1] = 0;         B[0][i*3+2] = 0;
            B[1][i*3+0] = 0;         B[1][i*3+1] = GradS[1];  B[1][i*3+2] = 0;
            B[2][i*3+0] = 0;         B[2][i*3+1] = 0;         B[2][i*3+2] = GradS[2];
            B[3][i*3+0] = GradS[1];  B[3][i*3+1] = GradS[0];  B[3][i*3+2] = 0;
            B[4][i*3+0] = 0;         B[4][i*3+1] = GradS[2];  B[4][i*3+2] = GradS[1];
            B[5][i*3+0] = GradS[2];  B[5][i*3+1] = 0;         B[5][i*3+2] = GradS[0];
        }

        return B;
    }

    inline helper::fixed_array<Deriv, 4> getShapeDerivatives(const Coord & a, const Coord & b, const Coord & c, const Coord & d) const
    {
        Real volume = fabs( dot( cross(b-a, c-a), d-a ) ) / 6.0;
        return {
            -(cross(b,c) + cross(c,d) + cross(d,b)) / (6*volume),
             (cross(c,d) + cross(d,a) + cross(a,c)) / (6*volume),
            -(cross(d,a) + cross(a,b) + cross(b,d)) / (6*volume),
             (cross(a,b) + cross(b,c) + cross(c,a)) / (6*volume)
        };
    }

    void draw(const core::visual::VisualParams* vparams) override;

protected:
    void computeStiffnessMatrix();

    inline Corotational corotational_method() const {
        helper::ReadAccessor<Data<helper::OptionsGroup>> m_corotational = d_corotational;
        unsigned int id = m_corotational->getSelectedId();
        switch (id) {
            case 0: return Corotational::NONE;
            case 1: return Corotational::LARGE;
            case 2: return Corotational::POLAR;
            case 3: return Corotational::SVD;
            case 4: return Corotational::DEFORMATION_GRADIENT;
            default:return Corotational::NONE;
        }
    }

    inline Mat33 extractRotationSVD(const Mat33 & A, const Mat33 & transformation, const Mat33 & initial_rotation) const {
        Mat33 F = A*transformation;

        if(determinant(F) < 1e-6) {
            // inverted or too flat element -> SVD decomposition + handle degenerated cases
            defaulttype::Mat<3,3,Real> U(defaulttype::NOINIT), V(defaulttype::NOINIT);
            defaulttype::Vec<3,Real> Sdiag(defaulttype::NOINIT);

            if (helper::Decompose<Real>::SVD_stable( F, U, Sdiag, V ))
                msg_warning() << "Rotation extraction of a degenerate element";

            return U * V.transposed() * initial_rotation;
        } else {
            // not inverted & not degenerated -> classical polar
            return extractRotationPolar(A);
        }
    }

    inline Mat33 extractRotationPolar(const Mat33 & A) const {
        Mat33 R(defaulttype::NOINIT);
        helper::Decompose<Real>::polarDecomposition( A, R );
        return R;
    }

    inline Mat33 extractRotationLarge(const Mat33 & A) const {
        const Coord & a = A[0]; // Edge 0
        const Coord & b = A[1]; // Edge 1

        Coord x = a.normalized();
        Coord y = b.normalized();
        Coord z = cross(x, y).normalized();

        return Mat33(x, y, z);
    }

    Mat33 computeRotation(const helper::fixed_array<Coord, 4> & node_position, const helper::fixed_array<Coord, 4> & node_rest_position, const TetrahedronID & tId) const;

    // Inputs
    Data< Real > d_youngModulus;
    Data< Real > d_poissonRatio;
    Data<Real> d_draw_scale; ///< Scale applied to tetrahedrons drawn when the draw forcefield option is enabled
    Data<helper::vector<Coord>> d_initial_positions;  ///< List of initial coordinates of the tetrahedrons nodes
    Data<helper::vector<Tetrahedron>> d_tetrahedrons; ///< List of tetrahedrons by their nodes indices (ex: [t1p1 t1p2 t1p3 t1p4 t2p1 t2p2 t2p3 t2p4...])
    Data<sofa::helper::OptionsGroup> d_corotational; ///< the computation method of the rotation extraction
    Data<bool> d_incremental_rotation; ///< Extract the incremental rotation (rotate the displacement with the last extracted rotation before the current extraction)

    // Outputs
    Data<helper::vector<Mat1212>> d_element_stiffness_matrices; ///< List of elements stiffness matrices
    Data<helper::vector<helper::fixed_array<Deriv, 4>>> d_element_shape_derivatives; ///< List of elements shape derivatives at their nodes

private:
    helper::vector<Mat33> m_element_initial_inverted_transformations; ///< Initial inverted transformation matrices extracted from each elements
    helper::vector<Mat33> m_element_initial_rotations; ///< Initial rotation matrices extracted from each elements
    helper::vector<Mat33> m_element_rotations; ///< Rotation matrices extracted from each elements
};

} // namespace forcefield

} // namespace caribou

} // namespace sofa

#endif //CARIBOU_FORCEFIELD_FEMFORCEFIELD_H
