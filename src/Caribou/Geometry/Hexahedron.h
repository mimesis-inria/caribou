#ifndef CARIBOU_GEOMETRY_HEXAHEDRON_H
#define CARIBOU_GEOMETRY_HEXAHEDRON_H

#include <Caribou/config.h>
#include <Caribou/Algebra/Vector.h>
#include <Caribou/Geometry/Node.h>
#include <Caribou/Geometry/Interpolation/Hexahedron.h>

namespace caribou {
namespace geometry {

template <typename CanonicalElementType>
struct Hexahedron : public CanonicalElementType
{
    static constexpr INTEGER_TYPE NumberOfNodes = CanonicalElementType::NumberOfNodes;
    using NodeType = caribou::geometry::Node<3>;
    using Index = std::size_t ;
    using Real = FLOATING_POINT_TYPE;

    constexpr
    Hexahedron()
            : p_nodes(CanonicalElementType::nodes)
    {}

    constexpr
    Hexahedron(const std::array<NodeType, NumberOfNodes> & nodes)
            : p_nodes(nodes)
    {}

    template <
            typename ...Nodes,
            REQUIRES(NumberOfNodes == sizeof...(Nodes)),
            REQUIRES(std::conjunction_v<std::is_same<NodeType, Nodes>...>)
    >
    constexpr
    Hexahedron(Nodes&&...remaining_nodes)
    : p_nodes {std::forward<Nodes>(remaining_nodes)...}
    {}

    /** Get the Node at given index */
    const NodeType &
    node(Index index) const
    {
        return p_nodes[index];
    }

    /** Get the Node at given index */
    NodeType &
    node(Index index)
    {
        return p_nodes[index];
    }

    /** Compute the jacobian matrix evaluated at local position {u,v,w}
     * (see interpolation::CanonicalElement::Jacobian for more details).
     * */
    inline algebra::Matrix<3, 3, Real>
    jacobian (const Real & u, const Real & v, const Real & w) const
    {
        return CanonicalElementType::Jacobian({u,v, w}, p_nodes);
    }

    /**
     * Compute an integral approximation by gauss quadrature on the hexahedron of the given evaluation function.
     *
     * @example
     * \code{.cpp}
     * // Integrate the polynomial 1 + 2x + 2xy + 3*z on an hexahedron.
     * float result = LinearHexahedron(x1, x2, x3, x4, x5, x6, x7, x8).gauss_integrate(
     *   [] (const LinearHexahedron & hexa, const float & xi, const float & eta, const float & zeta) -> float {
     *     return 1 + 2*xi + 2*xi*eta + 3*zeta;
     *   }
     * );
     * \endcode
     *
     * @tparam ValueType The result type of the evaluation function.
     * This type must implement the assignment "=", assignment-addition "+=", and multiplication "*" with a scalar type (float, double) operators.
     * @tparam EvaluateFunctionType Callback function reference type. See evaluate parameter.
     *
     * @param f
     * Callback function of the signature
     *
     *     ValueType f (const LinearHexahedron & hexa, const float & u, const float & v, const float & w);
     *
     * Where hexa is a reference to the current hexahadron on which we integrate, and the coordinates u, v and w
     * forms the local position of a sample point on which we want to get the evaluation value of type ValueType.
     *
     * @return The value of the integral computed on this hexahedron.
     *
     */
    template <typename ValueType , typename EvaluateFunctor>
    inline
    ValueType gauss_quadrature(const ValueType & initial_value, EvaluateFunctor f) const
    {
        static_assert(CanonicalElementType::gauss_nodes.size() == CanonicalElementType::gauss_weights.size(),
                "Gauss nodes must have assigned weights.");

        ValueType result = initial_value;

        for (std::size_t i = 0; i < CanonicalElementType::gauss_nodes.size(); ++i) {
            const auto p = CanonicalElementType::gauss_nodes[i];
            const auto w = CanonicalElementType::gauss_weights[i];
            const auto detJ = jacobian(p[0], p[1], p[2]).determinant();
            const auto eval = f(*this, p[0], p[1], p[2]);
            result += eval * w * detJ;
        }

        return result;
    }

    /**
     * Compute an integral approximation by gauss quadrature on the hexahedron of the given evaluation function.
     *
     * @example
     * \code{.cpp}
     * // Integrate the polynomial 1 + 2x + 2xy + 3*z on an hexahedron.
     * float result = Hexahedron(x1, x2, x3, x4, x5, x6, x7, x8).gauss_integrate(
     *   [] (const LinearHexahedron & hexa, const float & xi, const float & eta, const float & zeta) -> float {
     *     return 1 + 2*xi + 2*xi*eta + 3*zeta;
     *   }
     * );
     * \endcode
     *
     * @tparam EvaluateFunctionType Callback function reference type. See f parameter.
     *
     * @param f
     * Callback function of the signature
     *
     *     ValueType f (const Hexahedron & hexa, const float & u, const float & v, const float & w);
     *
     * Where hexa is a reference to the current hexahadron on which we integrate, and the coordinates u, v and w
     * forms the local position of a sample point on which we want to get the evaluation value.
     *
     * @return The value of the integral computed on this hexahedron.
     *
     */
    template <typename EvaluateFunctor>
    inline
    auto gauss_quadrature(EvaluateFunctor f) const
    {
        static_assert(CanonicalElementType::gauss_nodes.size() == CanonicalElementType::gauss_weights.size(),
                      "Gauss nodes must have assigned weights.");

        const auto p0 = CanonicalElementType::gauss_nodes[0];
        const auto w0 = CanonicalElementType::gauss_weights[0];
        const auto detJ0 = jacobian(p0[0], p0[1], p0[2]).determinant();
        const auto eval0 = f(*this, p0[0], p0[1], p0[2]);
        auto result = eval0 * w0 * detJ0;

        for (std::size_t i = 1; i < CanonicalElementType::gauss_nodes.size(); ++i) {
            const auto p = CanonicalElementType::gauss_nodes[i];
            const auto w = CanonicalElementType::gauss_weights[i];
            const auto detJ = jacobian(p[0], p[1], p[2]).determinant();
            const auto eval = f(*this, p[0], p[1], p[2]);
            result += eval * w * detJ;
        }

        return result;
    }

    /**
     * Check whether the hexahedron is a parallelepiped.
     *
     * A parallelepiped:
     * - Has six faces, each of which is a parallelogram
     * - Has three pairs of parallel faces
     *
     * The transformation from its elemental frame to its world frame can be defined as
     *                      | x1 + 1/2 (1 + u) hx |
     * (x,y,z) = Q(u,v,w) = | y1 + 1/2 (1 + v) hy |
     *                      | z1 + 1/2 (1 + w) hz |
     *
     * where (x1,y1,z1) are the world coordinates of the node #0 on the hexahedron, and (hx,hy,hz) denote the hexahedron
     * size w.r.t the x,y and z directions.
     *
     * And the Jacobian of this transformation is constant and defined as
     *                          | hx 0  0  |
     * J = gradQ^T(u,v,w) = 1/2 | 0  hy 0  |
     *                          | 0  0  hz |
     */
    inline
    bool
    is_a_parallelepiped() const noexcept
    {
        const auto lx = node(1) - node(0);
        const auto lz = node(4) - node(0);

        return (
                (node(3)-node(2) + lx).length_squared() < EPSILON && // edges 0-1 and 2-3 are of the same length
                (node(1)-node(5) + lz).length_squared() < EPSILON && // edges 0-4 and 1-5 are of the same length
                (node(2)-node(6) + lz).length_squared() < EPSILON && // edges 0-4 and 2-6 are of the same length
                (node(3)-node(7) + lz).length_squared() < EPSILON    // edges 0-4 and 3-7 are of the same length
        );
    }

private:
    std::array<NodeType, NumberOfNodes> p_nodes;
};

template <
        typename ...Nodes,
        REQUIRES(8 == sizeof...(Nodes)),
        REQUIRES(std::conjunction_v<std::is_same<caribou::geometry::Node<3>, Nodes>...>)
>
Hexahedron(Nodes&&...remaining_nodes) -> Hexahedron<interpolation::Hexahedron8>;

Hexahedron() -> Hexahedron<interpolation::Hexahedron8>;

} // namespace geometry
} // namespace caribou
#endif //CARIBOU_GEOMETRY_HEXAHEDRON_H
