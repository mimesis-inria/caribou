#ifndef CARIBOU_MECHANICS_ELASTICITY_STRAIN_H
#define CARIBOU_MECHANICS_ELASTICITY_STRAIN_H

#include <Caribou/config.h>
#include <Caribou/Traits.h>
#include <Eigen/Core>

namespace caribou::mechanics::elasticity::strain {

using Float = FLOATING_POINT_TYPE;

template<int nRows, int nColumns, int Options=0>
using Matrix = Eigen::Matrix<FLOATING_POINT_TYPE, nRows, nColumns, Options>;

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
template <int NumberOfNodes, int Dimension, int Options>
static inline
Matrix<Dimension, Dimension>
F (const Matrix<NumberOfNodes, Dimension, Options> & dN_dx, const Matrix<NumberOfNodes, Dimension, Options> & U)
{
    const auto Id = Matrix<Dimension, Dimension>::Identity();
    Matrix<Dimension, Dimension> GradU = U.row(0).transpose() * dN_dx.row(0);
    for (size_t i = 1; i < NumberOfNodes; ++i) {
        GradU.noalias() += U.row(i).transpose() * dN_dx.row(i);
    }
    return GradU + Id;
}

} // namespace caribou::mechanics::elasticity::strain

#endif //CARIBOU_MECHANICS_ELASTICITY_STRAIN_H
