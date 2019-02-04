#ifndef CARIBOU_MECHANICS_ELASTICITY_STRAIN_H
#define CARIBOU_MECHANICS_ELASTICITY_STRAIN_H

#include <Caribou/config.h>
#include <Caribou/Traits.h>
#include <Caribou/Algebra/Matrix.h>
#include <Caribou/Algebra/Vector.h>

namespace caribou {
namespace mechanics {
namespace elasticity {
namespace strain {

using Float = FLOATING_POINT_TYPE;

using Mat33 = algebra::Matrix<3,3, FLOATING_POINT_TYPE>;
using Vec3  = algebra::Matrix<3,1, FLOATING_POINT_TYPE>;

/**
 * Strain-displacement matrix B evaluated at local coordinates {xi, eta, zeta}.
 *
 * This function computes the strain-displacement matrix which regroups the partial derivatives of the shape functions
 * with respect to the world coordinates.
 *
 *       | dN0_dx    0       0          dN1_dx    0       0                    dNn_dx    0       0    |
 *       |   0     dN0_dy    0            0     dN1_dy    0                      0     dNn_dy    0    |
 * B   = |   0       0     dN0_dz         0       0     dN1_dz       ...         0       0     dNn_dz |
 *       | dN0_dy  dN0_dx    0          dN1_dy  dN1_dx    0                    dNn_dy  dNn_dx    0    |
 *       |   0     dN0_dz  dN0_dy         0     dN1_dz  dN1_dy                   0     dNn_dz  dNn_dy |
 *       | dN0_dz    0     dN0_dx       dN1_dz    0     dN1_dx                 dNn_dz    0     dNn_dx |
 *
 * with dNi_dx the ith derivatives of the shape function w.r.t the world coordinate x, and n the number of nodes (shape
 * functions) in the element.
 *
* @tparam ElementType
 *    The class type of the element on which we are computing the deformation gradient. The class must have the
 *    following public function:
 *       - jacobian(xi, eta, zeta)  -> Mat33 : Computes the jacobian matrix at local coordinates {xi, eta, zeta}
 *       - dN(xi, eta, zeta) -> Matrix<n, 3>  : Computes the derivatives of the n shape functions evaluated at local
 *         coordinates {xi, eta, zeta}
 *       - node(i) -> Vec3 : Get the world (global) coordinates of the node i
 *
 *    and the following public static member:
 *       - NumberOfNodes : The number of nodes this element type contains
 */
template <typename ElementType, typename LocalCoordinates, REQUIRES(ElementType::CanonicalDimension == 3)>
static inline
algebra::Matrix<6, ElementType::NumberOfNodes*3, FLOATING_POINT_TYPE>
B (const ElementType & element, const LocalCoordinates & coordinates)
{
    const auto J = element.jacobian(coordinates); // Jacobian matrix of the shape function N
    const auto Jinv = J^-1;
    const auto dN = ElementType::dN(coordinates);

    algebra::Matrix<6, ElementType::NumberOfNodes*3, FLOATING_POINT_TYPE> r; // Result (B)

    for (size_t i = 0; i < ElementType::NumberOfNodes; ++i) {
        const auto dN_dXi = dN.row(i).T();
        const auto dN_dx = Jinv.T() * dN_dXi;

        r(0, i*3+0) = dN_dx[0];  r(0, i*3+1) = 0;         r(0, i*3+2) = 0;
        r(1, i*3+0) = 0;         r(1, i*3+1) = dN_dx[1];  r(1, i*3+2) = 0;
        r(2, i*3+0) = 0;         r(2, i*3+1) = 0;         r(2, i*3+2) = dN_dx[2];
        r(3, i*3+0) = dN_dx[1];  r(3, i*3+1) = dN_dx[0];  r(3, i*3+2) = 0;
        r(4, i*3+0) = 0;         r(4, i*3+1) = dN_dx[2];  r(4, i*3+2) = dN_dx[1];
        r(5, i*3+0) = dN_dx[2];  r(5, i*3+1) = 0;         r(5, i*3+2) = dN_dx[0];
    }

    return r;
}

/**
 * Deformation gradient tensor F.
 *
 * @tparam NumberOfNodes Number of nodes inside the element
 * @tparam Dimension Dimension of every positions
 * @tparam DataType Floating point data type used
 *
 * @param dN_dx Matrix NxD where the ith row contains the ith derivative of the shape function w.r.t the jth world
 * coordinate unit (jth column).
 *
 * Example in 3D:
 *
 *         |dN0_dx dN0_dy dN0_dz|
 * dN_dx = |dN1_dx dN1_dy dN1_dz|
 *         |dN2_dx dN2_dy dN2_dz|
 *         |       ...          |
 *         |dNn_dx dNn_dy dNn_dz|
 *
 * @param U Matrix where the ith row contains the displacement vector of the ith element's node
 *
 *  Example in 3D:
 *
 *     |(x0 - X0)   (y0 - Y0)   (z0 - Z0)|
 *     |(x1 - X1)   (y1 - Y1)   (z1 - Z1)|
 * U = |(x2 - X2)   (y2 - Y2)   (z2 - Z2)|
 *     |                ...              |
 *     |(xn - Xn)   (yn - Yn)   (zn - Zn)|
 *
 * @return  DxD deformation gradient tensor F where D is the dimension of the world coordinates.
 */
template <size_t NumberOfNodes, size_t Dimension, typename DataType>
static inline
algebra::Matrix<Dimension, Dimension, DataType>
F (const algebra::Matrix<NumberOfNodes, Dimension, DataType> & dN_dx, const algebra::Matrix<NumberOfNodes, Dimension, DataType> & U)
{
    const Mat33 I = Mat33::Identity();
    Mat33 GradU = dN_dx.row(0).T() * U.row(0);
    for (size_t i = 1; i < NumberOfNodes; ++i) {
        GradU += dN_dx.row(i).T() * U.row(i);
    }
    return GradU + I;
}

/**
 * Deformation gradient tensor F evaluated at local coordinates {xi, eta, zeta}.
 *
 * This function computes the deformation gradient tensor which regroups the partial derivative of the displacement
 * vector with respect to the material coordinates (initial or undeformed configuration).
 *
 *     | du/dx  du/dy  dy/dz  |   | 1 0 0  |
 * F = | dv/dx  dv/dy  dv/dz  | + | 0 1 0  |
 *     | dw/dx  dw/dy  dw/dz  |   | 0 0 1  |
 *
 * with the displacement u (xi, eta, zeta) = {u, v, w} which is approximated from the node displacements of the element.
 *
 * @tparam ElementType
 *    The class type of the element on which we are computing the deformation gradient. The class must have the
 *    following public function:
 *       - jacobian(xi, eta, zeta)  -> Mat33 : Computes the jacobian matrix at local coordinates {xi, eta, zeta}
 *       - dN(xi, eta, zeta) -> Matrix<n, 3>  : Computes the derivatives of the n shape functions evaluated at local
 *         coordinates {xi, eta, zeta}
 *       - node(i) -> Vec3 : Get the world (global) coordinates of the node i
 *
 *    and the following public static member:
 *       - NumberOfNodes : The number of nodes this element type contains
 *
 */
template <typename ElementType, typename LocalCoordinates>
static inline Mat33
F (const ElementType & initial_element, const ElementType & deformed_element,
   const LocalCoordinates & coordinates)
{
    const Mat33 J = initial_element.jacobian(coordinates); // Jacobian matrix of the shape function N
    const Mat33 Jinv = J^-1;
    const Mat33 I = Mat33::Identity();

    Mat33 GradU;
    GradU.fill(0);
    const auto dN = ElementType::dN(coordinates);
    for (size_t i = 0; i < ElementType::NumberOfNodes; ++i) {
        const auto dN_dXi = dN.row(i).T();
        const auto dN_dx = Jinv.T() * dN_dXi;
        const auto u = deformed_element.node(i) - initial_element.node(i);
        GradU += dN_dx * u.T();
    }

    return GradU + I;
}

/**
 * Small (infinitesimal) strain tensor evaluated at local coordinates {xi, eta, zeta}.
 *
 * This function computes the small strain tensor epsilon :
 *
 * epsilon = 1/2 (F^T + F) - I
 *
 * where F is the deformation gradient tensor.
 *
 * @tparam ElementType
 *    see caribou::mechanics::elasticity::strain::F::ElementType
 */
template <typename ElementType, typename LocalCoordinates>
static inline Mat33
small_strain (const ElementType & initial_element, const ElementType & deformed_element,
              const LocalCoordinates & coordinates)
{
    const auto F = strain::F(initial_element, deformed_element, coordinates);
    const auto Ft = F.T();
    const auto I = Mat33::Identity();

    return 1./2 * (Ft + F) - I;
}

/**
 * Lagrangian finite strain tensor (also called the Green-Lagrangian strain tensor or Green – St-Venant strain tensor)
 * evaluated at local coordinates {xi, eta, zeta}.
 *
 * This function computes the strain tensor E :
 *
 * E = 1/2 (C - I)
 *
 * C = F^T * F
 *
 * where C is the right Cauchy–Green deformation tensor and F is the deformation gradient tensor.
 *
 * @tparam ElementType
 *    see caribou::mechanics::elasticity::strain::F::ElementType
 */
template <typename ElementType, typename LocalCoordinates>
static inline Mat33
strain (const ElementType & initial_element, const ElementType & deformed_element,
        const LocalCoordinates & coordinates)
{
    const auto F = strain::F(initial_element, deformed_element, coordinates);
    const auto Ft = F.T();
    const auto I = Mat33::Identity();
    const auto C = Ft*F; // right Cauchy–Green deformation tensor

    return 1./2 * (C - I);
}

} // namespace strain
} // namespace elasticity
} // namespace mechanics
} // namespace caribou

#endif //CARIBOU_MECHANICS_ELASTICITY_STRAIN_H
