#ifndef CARIBOU_GEOMETRY_TRIANGLE_H
#define CARIBOU_GEOMETRY_TRIANGLE_H

#include <Caribou/config.h>
#include <Caribou/Algebra/Vector.h>
#include <Caribou/Geometry/Node.h>
#include <Caribou/Geometry/Interpolation/Triangle.h>
#include <Caribou/Geometry/Internal/BaseTriangle.h>

namespace caribou {
namespace geometry {

template <size_t Dim, typename CanonicalElementType = interpolation::Triangle3>
struct Triangle : public internal::BaseTriangle<Dim, CanonicalElementType, Triangle<Dim, CanonicalElementType>>
{
    static constexpr INTEGER_TYPE NumberOfNodes = CanonicalElementType::NumberOfNodes;
    using NodeType = caribou::geometry::Node<Dim>;
    using Index = std::size_t ;
    using Real = FLOATING_POINT_TYPE;

    static_assert(Dim == 2 or Dim == 3, "Only 2D and 3D triangles are supported.");

    template <
            typename ...Nodes,
            REQUIRES(NumberOfNodes == sizeof...(Nodes)+1)
    >
    constexpr
    Triangle(const NodeType & first_node, Nodes&&...remaining_nodes)
    : p_nodes {first_node, std::forward<Nodes>(remaining_nodes)...}
    {}

    constexpr
    const NodeType &
    node(Index index) const
    {
        return p_nodes[index];
    }

    constexpr
    NodeType &
    node(Index index)
    {
        return p_nodes[index];
    }

    /** Compute the center position **/
    auto
    center() const noexcept
    {
        return (node(0) + node(1) + node(2)) / 3.0;
    }

    /** Compute the surface area **/
    FLOATING_POINT_TYPE
    area() const noexcept
    {
        auto n1 = node(0);
        auto n2 = node(1);
        auto n3 = node(2);

        if constexpr (Dim == 2) {
            caribou::algebra::Matrix<3,3,Real> M ({
                    {(Real) n1[0], (Real) n2[0], (Real) n3[0]},
                    {(Real) n1[1], (Real) n2[1], (Real) n3[1]},
                    {(Real)  1.,   (Real)  1.,   (Real)  1.}
            });

            return 1 / 2. * std::abs(M.determinant());
        } else {
            auto v1 = n3 - n1;
            auto v2 = n2 - n1;

            return v1.cross(v2).length() / 2.;
        }
    }

    /** Compute the jacobian matrix evaluated at local position {u,v} */
    algebra::Matrix<Dim, 2, Real>
    jacobian (const Real & u, const Real & v) const
    {
        return CanonicalElementType::Jacobian({u,v}, p_nodes);
    }

private:
    std::array<NodeType, NumberOfNodes> p_nodes;
};

} // namespace geometry
} // namespace caribou
#endif //CARIBOU_GEOMETRY_TRIANGLE_H
