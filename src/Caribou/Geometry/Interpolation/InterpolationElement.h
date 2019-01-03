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
 * such as a linear Quad). All the positions used in this element are specified by local coordinates (such as
 * barycentric coordinates) relative to the element frame axis.
 *
 * @tparam Dim The dimension of the interpolation element (this is the dimension of the local frame. If the element is a
 * 3D triangle, the dimension is still 2D as the local frame of the element is in 2D.
 * @tparam NNodes The number of interpolation nodes of the element. Usually, a linear element will only have an
 * interpolation node per corners. Higher degree elements will usually have additional nodes between its corners.
 * @tparam ElementType_ The explicit element type (that will inherit this base class).
 */
template<INTEGER_TYPE Dim, INTEGER_TYPE NNodes, typename ElementType_>
struct InterpolationElement
{
    using ElementType = ElementType_;
    static constexpr INTEGER_TYPE Dimension = Dim;
    static constexpr INTEGER_TYPE NumberOfNodes = NNodes;

    /**
     * Compute the Jacobian matrix evaluated at local coordinates.
     *
     * The Jacobian is defined as:
     *
     * In 1D:
     *
     * J(u) = dx/du = sum dNi/du * u
     *
     * In 2D:
     *
     *          | dx/dXi    dy/dXi  |   | sum dNi/dXi   u    sum dNi/dXi   v |
     * J(u,v) = | dx/dEta   dy/dEta | = | sum dNi/dEta  u    sum dNi/dEta  v |
     *
     * In 3D:
     *
     *            | dx/dXi    dy/dXi   dz/dXi   |   | sum dNi/dXi   u    sum dNi/dXi   v    sum dNi/dXi   w |
     * J(u,v,w) = | dx/dEta   dy/dEta  dz/dEta  | = | sum dNi/dEta  u    sum dNi/dEta  v    sum dNi/dEta  w |
     *            | dx/dZeta  dy/dZeta dz/dZeta |   | sum dNi/dZeta u    sum dNi/dZeta v    sum dNi/dZeta w |
     *
     * where dNi/dXi (resp. dEta and dZeta) is the partial derivative of the shape function at node i (ui, vi, wi)
     * w.r.t the local frame Xi {Xi, Eta, Zeta} evaluated at local coordinate  {u, v, w}.
     *
     * @example
     * \code{.cpp}
     * // Computes the 1D Jacobian evaluated at local coordinates 0.5
     * double J = LagrangeElement::Jacobian (0.5);
     *
     * // Computes the 2D Jacobian evaluated at local coordinates {-0.4, 0.2}
     * Matrix<2,2> J = LagrangeElement::Jacobian (-0.4, 0.2);
     *
     * // Computes the 3D Jacobian evaluated at local coordinates {-0.4, 0.2, 0.1}
     * Matrix<3,3> J = LagrangeElement::Jacobian (-0.4, 0.2, 0.1);
     * \endcode
     */
    template<typename ...Coordinates, REQUIRES(Dimension == sizeof...(Coordinates))>
    inline
    auto
    Jacobian (Coordinates &&...e) const
    {
        const auto local_coordinates = caribou::algebra::Vector<Dimension> {std::forward<Coordinates>(e)...};
        const auto shape_derivatives = self().dN(std::forward<Coordinates>(e)...);

        auto sum_of_derivatives = shape_derivatives.row(0);
        for (std::size_t i = 1; i < NumberOfNodes; ++i)
            sum_of_derivatives += shape_derivatives.row(i);

        if CONSTEXPR_IF (Dimension == 1) {
            return local_coordinates[0] * sum_of_derivatives[0];
        } else if CONSTEXPR_IF (Dimension == 2) {
            const auto D = caribou::algebra::Matrix {
                    sum_of_derivatives[0], sum_of_derivatives[0],
                    sum_of_derivatives[1], sum_of_derivatives[1]
            };
            return D.direct_multiplication(caribou::algebra::RowVector<2>{
                    local_coordinates[0], local_coordinates[1]
            });
        } else { // Dimension == 3
            const auto D = caribou::algebra::Matrix {
                    sum_of_derivatives[0], sum_of_derivatives[0], sum_of_derivatives[0],
                    sum_of_derivatives[1], sum_of_derivatives[1], sum_of_derivatives[1],
                    sum_of_derivatives[2], sum_of_derivatives[2], sum_of_derivatives[2]
            };
            return D.direct_multiplication(caribou::algebra::RowVector<3>{
                    local_coordinates[0], local_coordinates[1], local_coordinates[2]
            });
        }
    }

    /**
     * Get the shape values evaluated at local coordinates {u, v}.
     */
    inline
    auto
    N (const algebra::Vector<Dimension, FLOATING_POINT_TYPE> & coordinates) const
    {
        if CONSTEXPR_IF (Dimension == 1)
            return self().N(coordinates[0]);
        else if CONSTEXPR_IF (Dimension == 2)
            return self().N(coordinates[0], coordinates[1]);
        else if CONSTEXPR_IF (Dimension == 3)
            return self().N(coordinates[0], coordinates[1], coordinates[2]);
    }

    /**
     * Get the shape derivatives w.r.t the local frame {dN/du, dN/dv} evaluated at local coordinates {u, v}.
     */
    inline
    auto
    dN (const algebra::Vector<Dimension, FLOATING_POINT_TYPE> & coordinates) const
    {
        if CONSTEXPR_IF (Dimension == 1)
            return self().dN(coordinates[0]);
        else if CONSTEXPR_IF (Dimension == 2)
            return self().dN(coordinates[0], coordinates[1]);
        else if CONSTEXPR_IF (Dimension == 3)
            return self().dN(coordinates[0], coordinates[1], coordinates[2]);
    }

    /**
     * Interpolate a value at local coordinates from the given interpolation node values
     * @tparam ValueType Type of the value to interpolate.
     * This type must implement the multiplication operator with a floating point value (scalar) : ValueType * scalar.
     */
    template <typename ValueType>
    inline
    auto
    interpolate_at_local_position (algebra::Vector<Dimension, FLOATING_POINT_TYPE> && coordinates, const algebra::Vector<NumberOfNodes, ValueType> & values) const
    {
        const auto shapes = N(std::forward<algebra::Vector<Dimension, FLOATING_POINT_TYPE>>(coordinates));
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
    interpolate_at_local_position (algebra::Vector<Dimension, FLOATING_POINT_TYPE> && coordinates, ValueType && v0, Values &&... v) const
    {
        const auto shapes = N(std::forward<algebra::Vector<Dimension, FLOATING_POINT_TYPE>>(coordinates));
        const algebra::Vector<NumberOfNodes, ValueType> values {std::forward<ValueType>(v0), std::forward<Values>(v)...};
        return  shapes.dot(values);
    }

private:
    const ElementType &self () const
    {
        return static_cast<const ElementType &>(*this);
    }
};

} // namespace interpolation
} // namespace geometry
} // namespace caribou
#endif //CARIBOU_GEOMETRY_INTERPOLATION_INTERPOLATIONELEMENT_H
