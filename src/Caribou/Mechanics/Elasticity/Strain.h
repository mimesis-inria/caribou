#ifndef CARIBOU_MECHANICS_ELASTICITY_STRAIN_H
#define CARIBOU_MECHANICS_ELASTICITY_STRAIN_H

#include <Caribou/config.h>
#include <Caribou/Algebra/Matrix.h>

namespace caribou {
namespace mechanics {
namespace elasticity {
namespace strain {

using Float = FLOATING_POINT_TYPE;

using Mat33 = algebra::Matrix<3,3, FLOATING_POINT_TYPE>;
using Vec3 = algebra::Matrix<3,1, FLOATING_POINT_TYPE>;

/**
 * Deformation gradient tensor evaluated at local coordinates {xi, eta, zeta}.
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
 *       - Jacobian(xi, eta, zeta)  -> Mat33 : Computes the jacobian matrix at local coordinates {xi, eta, zeta}
 *       - dN_dXi(i, xi, eta, zeta) -> Vec3  : Computes the derivatives of the node i shape function evaluated at local
 *         coordinates {xi, eta, zeta}
 *       - node(i) -> Vec3 : Get the world (global) coordinates of the node i
 *
 *    and the following public static member:
 *       - NumberOfNodes : The number of nodes this element type contains
 *
 */
 template <typename ElementType>
inline Mat33
F (const ElementType & initial_element, const ElementType & deformed_element,
   const Float & xi, const Float & eta, const Float & zeta)
{
    const Mat33 J = initial_element.Jacobian(xi, eta, zeta); // Jacobian matrix of the shape function N
    const Mat33 Jinv = J^-1;
    const Mat33 I = Mat33::Identity();

    Mat33 GradU (true /* initialized_to_zero */);
    for (size_t i = 0; i < ElementType::NumberOfNodes; ++i) {
        const auto dN_dXi = initial_element.dN_dXi(i, xi, eta, zeta);
        const auto dN_dx = Vec3(Jinv * dN_dXi);
        const auto u = Vec3(deformed_element.node(i) - initial_element.node(i));
        GradU += dN_dx * u.transposed();
    }

    return GradU + I;
}

} // namespace strain
} // namespace elasticity
} // namespace mechanics
} // namespace caribou

#endif //CARIBOU_MECHANICS_ELASTICITY_STRAIN_H
