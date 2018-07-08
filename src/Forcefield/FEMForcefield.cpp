#include <sofa/core/ObjectFactory.h>
#include "FEMForcefield.h"
#include "../Helper/LinearAlgebra.h"

namespace sofa
{

using namespace core::objectmodel;
using namespace core::behavior;
using namespace component::topology;
using namespace defaulttype;

namespace caribou
{

namespace forcefield
{

template<class DataTypes>
FEMForcefield<DataTypes>::FEMForcefield()
    : d_youngModulus(initData(&d_youngModulus, Real(1000), "youngModulus", "Young's modulus of the material"))
    , d_poissonRatio(initData(&d_poissonRatio, Real(0.3),  "poissonRatio", "Poisson's ratio of the material"))
    , d_initial_positions(initData(&d_initial_positions,   "initial_positions", "List of initial coordinates of the tetrahedrons nodes"))
    , d_tetrahedrons(initData(&d_tetrahedrons, "tetrahedrons", "List of tetrahedrons by their nodes indices (ex: [t1p1 t1p2 t1p3 t1p4 t2p1 t2p2 t2p3 t2p4...])"))
    , d_corotational(initData(&d_corotational, "corotational","Corotation \"small\", \"large\" (by QR), \"polar\" or \"svd\" displacements"))

    , d_use_centroid_deformation_for_rotation_extraction(
        initData(&d_use_centroid_deformation_for_rotation_extraction, false, "use_centroid_deformation_for_rotation_extraction",
                 "Use the centroid point of an element to extract the rotation motion"))
{
    sofa::helper::WriteAccessor<Data<sofa::helper::OptionsGroup>> corotationalOptions = d_corotational;
    corotationalOptions->setNames(4,"NONE", "LARGE", "POLAR", "SVD"); // .. add other options
    corotationalOptions->setSelectedItem(0);
}

template<class DataTypes>
void FEMForcefield<DataTypes>::init()
{
    Inherit::init();

    // If no tetrahedrons indices set, get the ones from the container (automatically found in the context)
    if (d_tetrahedrons.getValue().empty()) {
        auto container = this->getContext()->template get<TetrahedronSetTopologyContainer>();
        if (container) {
            if (container->getNumberOfTetrahedra() == 0) {
                msg_error()
                    << "A tetrahedron topology container was found, but contained an empty set of tetrahedrons.";
            } else {
                d_tetrahedrons.setParent(&container->getTetrahedronDataArray());
            }
        } else {
            msg_error()
                << "A set of tetrahedrons or a tetrahedron topology containing a set of tetrahedrons must be provided.";
        }
    }

    // If no rest positions, get the ones from the mechanical object (automatically found in the context)
    if (d_initial_positions.getValue().empty()) {
        MechanicalState<DataTypes> * state = this->getContext()->template get<MechanicalState<DataTypes>>();
        if (state) {
            if (state->getSize() == 0) {
                msg_error() << "A mechanical object was found, but contained no initial coordinates.";
            } else {
                if (state->findData("rest_position")) {
                    d_initial_positions.setParent(state->findData("rest_position"));
                } else {
                    sofa::helper::ReadAccessor<Data<VecCoord>> X = state->readRestPositions();
                    sofa::helper::WriteAccessor<Data<VecCoord>> initial_positions = d_initial_positions;
                    initial_positions.resize(X.size());
                    for (PointID pid = 0; pid < X.size(); ++pid) {
                        initial_positions[pid] = X[pid];
                    }
                }
            }
        } else {
            msg_error()
                << "A set of initial coordinates must be provided, or a mechanical object must be found in the context.";
        }
    }

    // todo (jnbrunet2000@gmail.com): This should be called when the input data changes or is initialized
    computeStiffnessMatrix();
}

template<class DataTypes>
void FEMForcefield<DataTypes>::computeStiffnessMatrix()
{
    sofa::helper::ReadAccessor<Data<VecCoord>> X = d_initial_positions;
    sofa::helper::ReadAccessor<Data<sofa::helper::vector<Tetrahedron>>> tetrahedrons = d_tetrahedrons;
    sofa::helper::WriteAccessor<Data<sofa::helper::vector<Mat1212>>> stiffness_matrices = d_element_stiffness_matrices;

    // Make sure we got some rest coordinates and tetrahedrons
    if (X.empty() || tetrahedrons.empty())
        return;

    // Make sure that the set of positions matches the nodes id of the tetrahedrons
    for (const Tetrahedron & tetrahedron: tetrahedrons) {
        for (const PointID & id : tetrahedron) {
            if (id > X.size()-1) {
                msg_error() << "The set of rest positions does not match the set of tetrahedrons "
                               "(some nodes got an id greater than the number of rest positions provided).";
                return;
            }
        }
    }

    // Lame's coefficients and the material matrix
    const Real & youngModulus = this->d_youngModulus.getValue();
    const Real & poissonRatio = this->d_poissonRatio.getValue();

    const Real l = youngModulus*poissonRatio / ((1 + poissonRatio)*(1 - 2*poissonRatio));
    const Real m = youngModulus / (2 * (1 + poissonRatio));
    const Real a = l + 2*m;

    Mat66 C;
    C[0][0] = a; C[0][1] = l; C[0][2] = l; C[0][3] = 0; C[0][4] = 0; C[0][5] = 0;
    C[1][0] = l; C[1][1] = a; C[1][2] = l; C[1][3] = 0; C[1][4] = 0; C[1][5] = 0;
    C[2][0] = l; C[2][1] = l; C[2][2] = a; C[2][3] = 0; C[2][4] = 0; C[2][5] = 0;
    C[3][0] = 0; C[3][1] = 0; C[3][2] = 0; C[3][3] = m; C[3][4] = 0; C[3][5] = 0;
    C[4][0] = 0; C[4][1] = 0; C[4][2] = 0; C[4][3] = 0; C[4][4] = m; C[4][5] = 0;
    C[5][0] = 0; C[5][1] = 0; C[5][2] = 0; C[5][3] = 0; C[5][4] = 0; C[5][5] = m;

    // Compute the elemental stiffness matrices and initialize the rotation matrices
    stiffness_matrices.resize(tetrahedrons.size());
    m_element_rotations.resize(tetrahedrons.size(), Mat33::Identity());
    m_element_initial_rotations.resize(tetrahedrons.size(), Mat33::Identity());
    m_element_initial_inverted_transformations.resize(tetrahedrons.size(), Mat33::Identity());

    size_t i = 0;
    for (const Tetrahedron & tetrahedron: tetrahedrons) {
        Coord
            a = X[tetrahedron[0]],
            b = X[tetrahedron[1]],
            c = X[tetrahedron[2]],
            d = X[tetrahedron[3]];

        Mat612 B;

        if (corotational_method() == Corotational::NONE || d_use_centroid_deformation_for_rotation_extraction.getValue()) {
            B = getStrainDisplacement(a,b,c,d);
        } else {
            // Transformation matrix
            Mat33 A(b-a, // Edge 0
                    c-a, // Edge 1
                    d-a  // Edge 2
            );

            // For SVD, we need the initial inverted transformation matrix to detect degenerate cases
            if (corotational_method() == Corotational::SVD)
                m_element_initial_inverted_transformations[i] = A.inverted();

            // We extract the initial rotation matrix from the steady state
            if (corotational_method() == Corotational::LARGE)
                m_element_initial_rotations[i] = extractRotationLarge(A);
            else
                m_element_initial_rotations[i] = extractRotationPolar(A);

            // Compute the rotated strain-displacement matrix
            if (corotational_method() == Corotational::LARGE)
                B = getStrainDisplacement(a-a, b-a, c-a, d-a, m_element_initial_rotations[i]);
            else
                B = getStrainDisplacement(a,b,c,d, m_element_initial_rotations[i]);
        }

        // Compute the stiffness matrix
        Real volume = fabs( dot( cross(b-a, c-a), d-a ) ) / 6.0;
        stiffness_matrices[i] = volume*B.transposed()*C*B;

        ++i;
    }
}

template<class DataTypes>
void FEMForcefield<DataTypes>::addForce(const core::MechanicalParams* mparams, Data<VecDeriv>& d_f, const Data<VecCoord>& d_x, const Data<VecDeriv>& d_v)
{
    SOFA_UNUSED(mparams);
    SOFA_UNUSED(d_v);

    sofa::helper::WriteAccessor<Data<VecDeriv>> f = d_f;
    sofa::helper::ReadAccessor<Data<VecCoord>> x = d_x;
    sofa::helper::ReadAccessor<Data<VecCoord>> X = d_initial_positions;
    sofa::helper::ReadAccessor<Data<sofa::helper::vector<Tetrahedron>>> tetrahedrons = d_tetrahedrons;
    sofa::helper::ReadAccessor<Data<sofa::helper::vector<Mat1212>>> stiffness_matrices = d_element_stiffness_matrices;

    for (size_t i = 0; i < tetrahedrons.size(); ++i) {
        const auto & tetrahedron = tetrahedrons[i];

        PointID // Node indices
            a = tetrahedron[0],
            b = tetrahedron[1],
            c = tetrahedron[2],
            d = tetrahedron[3];

        Coord // Initial positions
            x0_a = X[a],
            x0_b = X[b],
            x0_c = X[c],
            x0_d = X[d];

        Coord // Current positions
            x_a = x[a],
            x_b = x[b],
            x_c = x[c],
            x_d = x[d];

        // Compute the current displacements for each nodes
        Deriv u[4];
        if (corotational_method() == Corotational::NONE) {
            // Linear method
            u[0]= (x_a - x0_a);
            u[1]= (x_b - x0_b);
            u[2]= (x_c - x0_c);
            u[3]= (x_d - x0_d);
        } else {
            // Linear corotational method
            Mat33 & R  = m_element_rotations[i];

            if (d_use_centroid_deformation_for_rotation_extraction.getValue()) {
                // Compute the centroid displacement gradient
                u[0]= (x_a - x0_a);
                u[1]= (x_b - x0_b);
                u[2]= (x_c - x0_c);
                u[3]= (x_d - x0_d);

                Real volume = fabs( dot( cross(x0_b-x0_a, x0_c-x0_a), x0_d-x0_a ) ) / 6.0;
                Deriv shape_derivatives[4] = {
                    -(cross(x0_b,x0_c) + cross(x0_c,x0_d) + cross(x0_d,x0_b)) / (6*volume),
                    (cross(x0_c,x0_d) + cross(x0_d,x0_a) + cross(x0_a,x0_c)) / (6*volume),
                    -(cross(x0_d,x0_a) + cross(x0_a,x0_b) + cross(x0_b,x0_d)) / (6*volume),
                    (cross(x0_a,x0_b) + cross(x0_b,x0_c) + cross(x0_c,x0_a)) / (6*volume)
                };

                Mat33 GradU;
                for (auto indice : tetrahedron) {
                    GradU += (u[indice] ^ shape_derivatives[indice]);
                }

                // Extract the rotation from the deformation gradient
                const Mat33 F = GradU + Mat33::Identity();
                R = extractRotationPolar(F);
                const Mat33 Rt = R.transposed();

                u[0]= (Rt*x_a - x0_a);
                u[1]= (Rt*x_b - x0_b);
                u[2]= (Rt*x_c - x0_c);
                u[3]= (Rt*x_d - x0_d);

            } else {
                // Transformation matrix
                Mat33 A(
                    R * (x_b - x_a),
                    R * (x_c - x_a),
                    R * (x_d - x_a)
                );

                // Extract the current rotation
                Mat33 Rtemp;

                if (corotational_method() == Corotational::SVD)
                    Rtemp = extractRotationSVD(A, m_element_initial_inverted_transformations[i],
                                               m_element_initial_rotations[i]);
                else if (corotational_method() == Corotational::LARGE)
                    Rtemp = extractRotationLarge(A);
                else
                    Rtemp = extractRotationPolar(A);

                R = Rtemp * R;

                // Get the initial rotation
                const Mat33 &R0 = m_element_initial_rotations[i];

                // Rotate the displacements
                if (corotational_method() == Corotational::LARGE) {
                    u[0] = Coord(0, 0, 0);

                    u[1] = (R * (x_b - x_a) - R0 * (x0_b - x0_a));
                    u[1][1] = 0;
                    u[1][2] = 0;

                    u[2] = (R * (x_c - x_a) - R0 * (x0_c - x0_a));
                    u[2][2] = 0;

                    u[3] = (R * (x_d - x_a) - R0 * (x0_d - x0_a));
                } else {
                    u[0] = (R * x_a - R0 * x0_a);
                    u[1] = (R * x_b - R0 * x0_b);
                    u[2] = (R * x_c - R0 * x0_c);
                    u[3] = (R * x_d - R0 * x0_d);
                }
            }
        }

        // Gather the nodal displacement into one vector
        const Vec12 U (
            u[0][0], u[0][1], u[0][2], u[1][0], u[1][1], u[1][2], u[2][0], u[2][1], u[2][2], u[3][0], u[3][1], u[3][2]
        );

        // Compute the internal elastic force
        const Mat1212 & K = stiffness_matrices[i];
        const Vec12 F = K*U;

        // Return the computed nodal forces
        if (corotational_method() == Corotational::NONE) {
            f[a] -= Deriv( F[0], F[1], F[2] );
            f[b] -= Deriv( F[3], F[4], F[5] );
            f[c] -= Deriv( F[6], F[7], F[8] );
            f[d] -= Deriv( F[9], F[10], F[11] );
        } else {
            if (d_use_centroid_deformation_for_rotation_extraction.getValue()) {
                const Mat33 R = m_element_rotations[i];
                f[a] -= R * Deriv(F[0], F[1], F[2]);
                f[b] -= R * Deriv(F[3], F[4], F[5]);
                f[c] -= R * Deriv(F[6], F[7], F[8]);
                f[d] -= R * Deriv(F[9], F[10], F[11]);
            } else {
                const Mat33 Rt = m_element_rotations[i].transposed();
                f[a] -= Rt * Deriv(F[0], F[1], F[2]);
                f[b] -= Rt * Deriv(F[3], F[4], F[5]);
                f[c] -= Rt * Deriv(F[6], F[7], F[8]);
                f[d] -= Rt * Deriv(F[9], F[10], F[11]);
            }
        }
    }

}

template<class DataTypes>
void FEMForcefield<DataTypes>::addDForce(const core::MechanicalParams* mparams, Data<VecDeriv>& d_df, const Data<VecDeriv>& d_dx)
{
    sofa::helper::WriteAccessor<Data<VecDeriv>> df = d_df;
    sofa::helper::ReadAccessor<Data<VecCoord>> dx = d_dx;
    sofa::helper::ReadAccessor<Data<sofa::helper::vector<Tetrahedron>>> tetrahedrons = d_tetrahedrons;
    sofa::helper::ReadAccessor<Data<sofa::helper::vector<Mat1212>>> stiffness_matrices = d_element_stiffness_matrices;
    Real kFactor = (Real) mparams->kFactorIncludingRayleighDamping(this->rayleighStiffness.getValue());

    for (size_t i = 0; i < tetrahedrons.size(); ++i) {
        const auto & indices = tetrahedrons[i];

        const Mat33 & R  = m_element_rotations[i];
        const Mat33   Rt = R.transposed();

        PointID a = indices[0];
        PointID b = indices[1];
        PointID c = indices[2];
        PointID d = indices[3];

        const Deriv u[4] = {
            Rt*dx[a],
            Rt*dx[b],
            Rt*dx[c],
            Rt*dx[d],
        };

        const Vec12 U (
            u[0][0],
            u[0][1],
            u[0][2],
            u[1][0],
            u[1][1],
            u[1][2],
            u[2][0],
            u[2][1],
            u[2][2],
            u[3][0],
            u[3][1],
            u[3][2]
        );

        const Mat1212 & K = stiffness_matrices[i];
        const Vec12 dF = -1. * (K*U) * kFactor;

        df[a] += R * Deriv( dF[0], dF[1], dF[2] );
        df[b] += R * Deriv( dF[3], dF[4], dF[5] );
        df[c] += R * Deriv( dF[6], dF[7], dF[8] );
        df[d] += R * Deriv( dF[9], dF[10], dF[11] );
    }
}

template<class DataTypes>
void FEMForcefield<DataTypes>::addKToMatrix(BaseMatrix * matrix, SReal kFact, unsigned int &offset)
{
    SOFA_UNUSED(offset);
    sofa::helper::ReadAccessor<Data<sofa::helper::vector<Tetrahedron>>> tetrahedrons = d_tetrahedrons;
    sofa::helper::ReadAccessor<Data<sofa::helper::vector<Mat1212>>> stiffness_matrices = d_element_stiffness_matrices;

    for (size_t i = 0; i < tetrahedrons.size(); ++i) {
        const auto &indices = tetrahedrons[i];
        const Mat1212 & stiffness = stiffness_matrices[i];

        for (size_t j = 0; j < indices.size(); ++j) {
            PointID idx_j = indices[j];
            for (size_t k = 0; k < indices.size(); ++k) {
                PointID idx_k = indices[k];

                Mat33 K;
                stiffness.getsub(j*3, k*3, K);
                if (corotational_method() == Corotational::NONE) {
                    K = -kFact * K;
                } else {
                    const Mat33 & R  = m_element_rotations[i];
                    const Mat33   Rt = R.transposed();

                    if (d_use_centroid_deformation_for_rotation_extraction.getValue()) {
                        K = -kFact * R * K * Rt;
                    } else {
                        K = -kFact * Rt * K * R;
                    }
                }

                for (unsigned l = 0; l < 3; ++l) {
                    for (unsigned m = 0; m < 3; ++m) {
                        matrix->add(idx_j*3 + l, idx_k*3 + m, K[l][m]);
                    }
                }
            }
        }
    }
}

SOFA_DECL_CLASS(FEMForcefield)
static int FEMForcefieldClass = core::RegisterObject("Caribou FEM Forcefield")
                                         .add< FEMForcefield<sofa::defaulttype::Vec3dTypes> >(true)
;


} // namespace forcefield

} // namespace caribou

} // namespace sofa