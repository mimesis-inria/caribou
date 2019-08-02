#ifndef CARIBOU_GEOMETRY_INTERPOLATION_CANONICALELEMENT_H
#define CARIBOU_GEOMETRY_INTERPOLATION_CANONICALELEMENT_H

#include <Caribou/config.h>
#include <Eigen/Core>

namespace caribou::geometry::interpolation {

/**
 * This is a base class that should be inherited by a explicit interpolation element types (example, a Lagrange element
 * such as a linear Quad). All the positions used in this element are specified by local coordinates relative to the
 * canonical element frame axis.
 *
 * @tparam Dim The dimension of the canonical element (this is the dimension of the local frame. If the element is a
 * 3D triangle, the dimension is still 2D since the canonical triangle is in 2D.
 * @tparam NNodes The number of interpolation nodes of the element. Usually, a linear element will only have an
 * interpolation node per corners. Higher degree elements will usually have additional nodes between its corners.
 * @tparam ElementType_ The explicit element type (that will inherit this base class).
 */
template<UNSIGNED_INTEGER_TYPE Dim, UNSIGNED_INTEGER_TYPE NNodes, typename CanonicalElementType_>
struct CanonicalElement
{
    static constexpr UNSIGNED_INTEGER_TYPE CanonicalDimension = Dim;
    static constexpr UNSIGNED_INTEGER_TYPE NumberOfNodes = NNodes;

    using CanonicalElementType = CanonicalElementType_;
    using LocalCoordinates = Eigen::Matrix<FLOATING_POINT_TYPE, 1, CanonicalDimension>;

    /**
     * Compute the Jacobian matrix evaluated at local coordinates.
     *
     * The Jacobian is defined as:
     *
     * 1D canonical element
     * --------------------
     *
     * 1D manifold:    J(u) = dx/du = sum_i dNi/du * xi
     *
     * 2D manifold:    J(u) = | dx/du | = | sum dNi/du xi |
     *                        | dy/du | = | sum dNi/du yi |
     *
     *                        | dx/du | = | sum dNi/du xi |
     * 3D manifold:    J(u) = | dy/du | = | sum dNi/du yi |
     *                        | dz/du | = | sum dNi/du zi |
     *
     *
     * 2D canonical element
     * --------------------
     *
     * 1D manifold:    N/A
     *
     * 2D manifold:    J(u,v) = | dx/du   dx/dv |   | sum dNi/du  xi    sum dNi/dv  xi |
     *                          | dy/du   dy/dv | = | sum dNi/du  yi    sum dNi/dv  yi |
     *
     *                          | dx/du   dx/dv |   | sum dNi/du  xi    sum dNi/dv  xi |
     * 3D manifold:    J(u,v) = | dy/du   dy/dv | = | sum dNi/du  yi    sum dNi/dv  yi |
     *                          | dz/du   dz/dv | = | sum dNi/du  zi    sum dNi/dv  zi |
     *
     *
     * 3D canonical element
     * --------------------
     *
     * 1D manifold:    N/A
     * 2D manifold:    N/A
     *                            | dx/du   dx/dv   dx/dw |   | sum dNi/du xi   sum dNi/dv xi    sum dNi/dw xi |
     * 3D manifold:    J(u,v,w) = | dy/du   dy/dv   dy/dw | = | sum dNi/du yi   sum dNi/dv yi    sum dNi/dw yi |
     *                            | dz/du   dz/dv   dz/dw |   | sum dNi/du zi   sum dNi/dv zi    sum dNi/dw zi |
     *
     *
     *
     *
     * where dNi/du (resp. dv and dw) is the partial derivative of the shape function at node i
     * w.r.t the local frame of the canonical element evaluated at local coordinate  {u, v, w} and
     * where {xi, yi and zi} are the world coordinates of the position of node i on its element manifold.
     *
     * @example
     * \code{.cpp}
     * // Computes the Jacobian of a 3D segment and its determinant evaluated at local coordinates 0.5 (half-way through the segment)
     * Segment<3> segment {{5, 5, 5}, {10, 5,0}};
     * Matrix<3,1> J = segment.Jacobian (0.5);
     * double detJ = (J^T * J).determinant() ^ 1/2;
     *
     * // Computes the Jacobian of a 3D triangle and its determinant evaluated at local coordinates {1/3, 1/3} (on its center point)
     * Triangle<3> triangle {{5,5,5}, {15, 5, 5}, {10, 10, 10}};
     * Matrix<3,2> J = triangle.Jacobian(1/3., 1/3.);
     * double detJ = (J^T * J).determinant() ^ 1/2;
     *
     * // Computes the Jacobian of a 3D tetrahedron and its determinant evaluated at local coordinates {1/3, 1/3, 1/3} (on its center point)
     * Tetrahedron<3> tetra {{5,5,5}, {15, 5, 5}, {10, 10, 10}};
     * Matrix<3,3> J = tetra.Jacobian(1/3., 1/3., 1/3.);
     * double detJ = J.determinant();
     * \endcode
     */
    template<typename LocalCoordinates, int WorldDimension>
    static inline
    Eigen::Matrix<FLOATING_POINT_TYPE, WorldDimension, CanonicalDimension>
    Jacobian (LocalCoordinates && coordinates, const Eigen::Matrix<FLOATING_POINT_TYPE, NumberOfNodes, WorldDimension> & nodes)
    {
        const auto shape_derivatives = CanonicalElementType::dL(std::forward<LocalCoordinates>(coordinates));

        auto Jt = Eigen::Matrix<FLOATING_POINT_TYPE, CanonicalDimension, WorldDimension>::Zero();
        auto positions = nodes.transpose();

        return positions * shape_derivatives;


//        if constexpr (CanonicalDimension == 1) { // Canonical element of dimension 1
//            auto sum = nodes[0] * shape_derivatives.row(0).transpose();
//            for (std::size_t i = 1; i < NumberOfNodes; ++i)
//                sum += nodes[i] * shape_derivatives.row(i).transpose();
//            return sum;
//        } else if constexpr (CanonicalDimension == 2) { // Canonical element of dimension 2
//
//            Eigen::Matrix<FLOATING_POINT_TYPE, WorldDimension, CanonicalDimension> J;
//
//
//
//            Matrix positions (
//                    nodes[0].T(), // [x, y, z]
//                    nodes[0].T()  // [x, y, z]
//            );
//            auto sum =
//                    positions
//                    .direct_multiplication(shape_derivatives.row(0).transposed()) // [du, dv]^T
//                    .transposed();
//            for (std::size_t i = 1; i < NumberOfNodes; ++i) {
//                positions = Matrix(
//                        nodes[i].T(), // [x, y, z]
//                        nodes[i].T()  // [x, y, z]
//                );
//                sum += positions.direct_multiplication(shape_derivatives.row(i).transposed()) // [du, dv]^T
//                        .transposed();
//            }
//            return sum;
//        } else { // Dimension == 3
//            Matrix positions (
//                    nodes[0].T(), // [x, y, z]
//                    nodes[0].T(), // [x, y, z]
//                    nodes[0].T()  // [x, y, z]
//            );
//            auto sum = positions
//                    .direct_multiplication(shape_derivatives.row(0).transposed()) // [du, dv, dw]^T
//                    .transposed();
//
//            for (std::size_t i = 1; i < NumberOfNodes; ++i) {
//                positions = Matrix(
//                        nodes[i].T(), // [x, y, z]
//                        nodes[i].T(), // [x, y, z]
//                        nodes[i].T()  // [x, y, z]
//                );
//                sum += positions.direct_multiplication(shape_derivatives.row(i).transposed()) // [du, dv, dw]^T
//                        .transposed();
//            }
//            return sum;
//        }
    }

    /**
     * Interpolate a value at local coordinates from the given interpolation node values specified as a matrix
     * @tparam ValueType Type of the value to interpolate.
     * @tparam ValueDimension Dimension of the value matrix.
     * This type must implement the multiplication operator with a floating point value (scalar) : ValueType * scalar.
     */
    template <typename ValueType, int ValueDimension, typename LocalCoordinates>
    static inline
    auto
    interpolate_at_local_position (LocalCoordinates && coordinates,
                                   const Eigen::Matrix<ValueType, NumberOfNodes, ValueDimension> & values)
    {
        return Eigen::Matrix<ValueType, ValueDimension, 1>(
            (values.array().colwise() * CanonicalElementType::L(std::forward<LocalCoordinates>(coordinates)).array()).matrix().colwise().sum().transpose()
            );
    }
};

} // namespace caribou::geometry::interpolation
#endif //CARIBOU_GEOMETRY_INTERPOLATION_CANONICALELEMENT_H
