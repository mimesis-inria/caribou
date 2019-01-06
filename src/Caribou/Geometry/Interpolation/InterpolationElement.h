#ifndef CARIBOU_GEOMETRY_INTERPOLATION_INTERPOLATIONELEMENT_H
#define CARIBOU_GEOMETRY_INTERPOLATION_INTERPOLATIONELEMENT_H

#include <Caribou/config.h>
#include <Caribou/Algebra/Matrix.h>
#include <Caribou/Algebra/Vector.h>

namespace caribou {
namespace geometry {
namespace interpolation {

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
template<typename ElementType_, typename CanonicalElementType_>
struct InterpolationElement
{
    using ElementType = ElementType_;
    using CanonicalElementType = CanonicalElementType_;
    static constexpr INTEGER_TYPE Dimension = CanonicalElementType::Dimension;
    static constexpr INTEGER_TYPE NumberOfNodes = CanonicalElementType::NumberOfNodes;

    using Real = FLOATING_POINT_TYPE;
    using LocalCoordinates = algebra::Vector<Dimension, FLOATING_POINT_TYPE>;

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
    template<typename ...Coordinates, REQUIRES(Dimension == sizeof...(Coordinates))>
    inline
    auto
    Jacobian (Coordinates &&...e) const
    {
        using namespace caribou::algebra;

        const auto shape_derivatives = dN(std::forward<Coordinates>(e)...);

        if CONSTEXPR_IF (Dimension == 1) { // Canonical element of dimension 1
            auto sum = self().node(0) * shape_derivatives.row(0).transposed();
            for (std::size_t i = 1; i < NumberOfNodes; ++i)
                sum += self().node(0) * shape_derivatives.row(i).transposed();
            return sum;
        } else if CONSTEXPR_IF (Dimension == 2) { // Canonical element of dimension 2
            Matrix positions (
                    self().node(0), // [x, y, z]
                    self().node(0)  // [x, y, z]
            );
            auto sum =
                    positions
                    .direct_multiplication(shape_derivatives.row(0).transposed()) // [du, dv]^T
                    .transposed();
            for (std::size_t i = 1; i < NumberOfNodes; ++i) {
                positions = Matrix(
                        self().node(i), // [x, y, z]
                        self().node(i)  // [x, y, z]
                );
                sum += positions.direct_multiplication(shape_derivatives.row(i).transposed()) // [du, dv]^T
                        .transposed();
            }
            return sum;
        } else { // Dimension == 3
            Matrix positions (
                    self().node(0), // [x, y, z]
                    self().node(0), // [x, y, z]
                    self().node(0)  // [x, y, z]
            );
            auto sum = positions
                    .direct_multiplication(shape_derivatives.row(0).transposed()) // [du, dv, dw]^T
                    .transposed();

            for (std::size_t i = 1; i < NumberOfNodes; ++i) {
                positions = Matrix(
                        self().node(i), // [x, y, z]
                        self().node(i), // [x, y, z]
                        self().node(i)  // [x, y, z]
                );
                sum += positions.direct_multiplication(shape_derivatives.row(i).transposed()) // [du, dv, dw]^T
                        .transposed();
            }
            return sum;
        }
    }

    /**
     * Get the shape values for each nodes evaluated at local coordinates {u, v}.
     */
    static constexpr
    algebra::Vector <NumberOfNodes, Real>
    N (const Real &u, const Real &v)
    {
        return get_N_shapes(u, v, std::make_index_sequence<NumberOfNodes>{});
    }

    /**
     * Get the shape values for each nodes evaluated at local coordinates {u, v}.
     */
    static constexpr
    algebra::Vector <NumberOfNodes, Real>
    N (const algebra::Vector<Dimension, FLOATING_POINT_TYPE> & coordinates)
    {
        if CONSTEXPR_IF (Dimension == 1)
            return N(coordinates[0]);
        else if CONSTEXPR_IF (Dimension == 2)
            return N(coordinates[0], coordinates[1]);
        else if CONSTEXPR_IF (Dimension == 3)
            return N(coordinates[0], coordinates[1], coordinates[2]);
    }

    /**
     * Get the shape derivatives for each nodes  w.r.t the local frame {dN/du, dN/dv} evaluated at local coordinates {u, v}.
     */
    static constexpr
    algebra::Matrix<NumberOfNodes, 2, Real>
    dN (const Real &u, const Real &v)
    {
        return get_N_shape_derivatives(u,v, std::make_index_sequence<NumberOfNodes>{});
    }


    /**
     * Get the shape derivatives for each nodes  w.r.t the local frame {dN/du, dN/dv} evaluated at local coordinates {u, v}.
     */
    static constexpr
    algebra::Matrix<NumberOfNodes, 2, Real>
    dN (const algebra::Vector<Dimension, FLOATING_POINT_TYPE> & coordinates)
    {
        if CONSTEXPR_IF (Dimension == 1)
            return dN(coordinates[0]);
        else if CONSTEXPR_IF (Dimension == 2)
            return dN(coordinates[0], coordinates[1]);
        else if CONSTEXPR_IF (Dimension == 3)
            return dN(coordinates[0], coordinates[1], coordinates[2]);
    }

    /**
     * Interpolate a value at local coordinates from the given interpolation node values
     * @tparam ValueType Type of the value to interpolate.
     * This type must implement the multiplication operator with a floating point value (scalar) : ValueType * scalar.
     */
    template <typename ValueType>
    inline
    auto
    interpolate_at_local_position (LocalCoordinates && position, const algebra::Vector<NumberOfNodes, ValueType> & values) const
    {
        const auto shapes = N(std::forward<LocalCoordinates>(position));
        return  shapes.dot(values);
    }

    /**
     * Interpolate a value at local coordinates from the given interpolation node values
     * @tparam ValueType Type of the value to interpolate.
     * This type must implement the multiplication operator with a floating point value (scalar) : ValueType * scalar.
     */
    template <typename ValueType, typename ...Values, REQUIRES(NumberOfNodes == sizeof...(Values)+1)>
    inline
    auto
    interpolate_at_local_position (LocalCoordinates && position, ValueType && v0, Values &&... v) const
    {
        const auto shapes = N(std::forward<LocalCoordinates>(position));
        const algebra::Vector<NumberOfNodes, ValueType> values {std::forward<ValueType>(v0), std::forward<Values>(v)...};
        return  shapes.dot(values);
    }

private:
    const ElementType &self () const
    {
        return static_cast<const ElementType &>(*this);
    }


    /**
     * Gather the N first shape functions (L) into a vector.
     * @example
     * \code{.cpp}
     * // Get the shape values at indices 0, 1, 2 and 3 into a vector.
     * Vector<4> shapes = LagrangeElement::get_N_shapes(u, v, 0, 1, 2, 3);
     * \endcode
     */
    template <std::size_t... Ix>
    static constexpr
    algebra::Vector <NumberOfNodes, Real>
    get_N_shapes (const Real &u, const Real &v, std::index_sequence<Ix...>)
    {
        return algebra::Vector <NumberOfNodes, Real> {
            CanonicalElementType::template L<Ix>(u, v)...
        };
    }

    /**
     * Gather the N first shape derivatives (dL) into a Matrix.
     * @example
     * \code{.cpp}
     * // Get the shape derivatives at indices 0, 1, 2 and 3 into a vector.
     * Matrix<4,2> shapes_derivatives = LagrangeElement::get_N_shape_derivatives(u, v, 0, 1, 2, 3);
     * \endcode
     */
    template <std::size_t... Ix>
    static constexpr
    algebra::Matrix<NumberOfNodes, 2, Real>
    get_N_shape_derivatives (const Real &u, const Real &v, std::index_sequence<Ix...>)
    {
        return algebra::Matrix<NumberOfNodes, 2, Real> {
            CanonicalElementType::template dL<Ix>(u, v)...
        };
    }
};

} // namespace interpolation
} // namespace geometry
} // namespace caribou
#endif //CARIBOU_GEOMETRY_INTERPOLATION_INTERPOLATIONELEMENT_H
