#ifndef CARIBOU_GEOMETRY_HEXAHEDRON_H
#define CARIBOU_GEOMETRY_HEXAHEDRON_H

#include <Eigen/Dense>

#include <Caribou/config.h>
#include <Caribou/Geometry/Quad.h>
#include <Caribou/Geometry/Interpolation/Hexahedron.h>
#include <Caribou/Geometry/Internal/BaseHexahedron.h>

namespace caribou::geometry {

template <typename CanonicalElementType>
struct Hexahedron : public internal::BaseHexahedron<CanonicalElementType, Hexahedron<CanonicalElementType>>
{
    static constexpr INTEGER_TYPE NumberOfNodes = CanonicalElementType::NumberOfNodes;

    using Base = internal::BaseHexahedron<CanonicalElementType, Hexahedron<CanonicalElementType>>;

    using LocalCoordinates = typename Base::LocalCoordinates;
    using WorldCoordinates = typename Base::WorldCoordinates;

    using BoundaryType = Quad<3, typename CanonicalElementType::BoundaryType>;

    template<int nRows, int nColumns, int Options=Eigen::RowMajor>
    using Matrix = Eigen::Matrix<FLOATING_POINT_TYPE, nRows, nColumns, Options>;

    template<int nRows, int nColumns>
    using Map = Eigen::Map<const Matrix<nRows, nColumns>>;

    template<int nRows, int Options=0>
    using Vector = Eigen::Matrix<FLOATING_POINT_TYPE, nRows, 1, Options>;

    template<int nRows>
    using MapVector = Eigen::Map<const Vector<nRows>>;

    Hexahedron()
    : p_nodes(Map<NumberOfNodes, 3>(&CanonicalElementType::nodes[0][0]))
    {
    }

    template <
            typename ...Nodes,
            REQUIRES(NumberOfNodes == sizeof...(Nodes)+1)
    >
    Hexahedron(const WorldCoordinates & first_node, Nodes&&...remaining_nodes)
    {
        construct_from_nodes<0>(first_node, std::forward<Nodes>(remaining_nodes)...);
    }

    Hexahedron(const Matrix<NumberOfNodes, 3> & m)
    : p_nodes(m)
    {}

    /** Get the Node at given index */
    inline
    auto
    node(UNSIGNED_INTEGER_TYPE index) const
    {
        return WorldCoordinates(p_nodes.row(index));
    }

    /** Get the Node at given index */
    inline
    auto
    node(UNSIGNED_INTEGER_TYPE index)
    {
        return WorldCoordinates(p_nodes.row(index));
    }

    /** Get a reference to the set of nodes */
    inline
    const auto &
    nodes() const
    {
        return p_nodes;
    }

    /**
     * Get the ith quadrangle face.
     */
    inline
    BoundaryType
    face(UNSIGNED_INTEGER_TYPE index) const
    {
        const auto & face_indices = CanonicalElementType::faces[index];

        Matrix<BoundaryType::NumberOfNodes, 3> m;
        for (std::size_t i = 0; i < BoundaryType::NumberOfNodes; ++i)
            m.row(i) = node(face_indices[i]);

        return QuadType(m);
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
                (node(3)-node(2) + lx).squaredNorm() < EPSILON && // edges 0-1 and 2-3 have the same length
                (node(1)-node(5) + lz).squaredNorm() < EPSILON && // edges 0-4 and 1-5 have the same length
                (node(2)-node(6) + lz).squaredNorm() < EPSILON && // edges 0-4 and 2-6 have the same length
                (node(3)-node(7) + lz).squaredNorm() < EPSILON    // edges 0-4 and 3-7 have the same length
        );
    }

    /**
     * Extract the frame positioned at the center of the hexahedron by computing the cross product of the unit vectors
     * from the center of the hexahedron to the center of the opposite faces.
     *
     * This function will return a matrix of the form:
     * | ux vx wx |
     * | uy vy wy |
     * | uz vz wz |
     *
     * Where [ux uy uz], [vx, vy, vz] and [wx, wy, wz] are orthogonal unitary vectors representing the u, v and w frame
     * in the current hexahedron. If the hexahedron is regular and not rotated, this matrix is the Identity matrix.
     * If it is regular but rotated, rotating the hexa by the transposed of this frame should align the u,v,w axis to the
     * x,y,z world frame (identity matrix).
     */
    inline
   Matrix<3, 3>
    frame(const LocalCoordinates & local_point) const
    {
        // Position of the point inside the hexahedron where the frame should be computed
        const auto p = T( local_point );

        // Project of the point on the quad facing the u axis
        const auto projected_on_u = T({1, local_point[1],local_point[2]});

        // Project of the point on the quad facing the v axis
        const auto projected_on_v = T({local_point[0], 1,local_point[2]});

        /* @todo(jnbrunet2000@gmail.com): select between the pairs of axis (center-to-u, center-to-v),
          (center-to-u, center-to-w) and (center-to-v, center-to-w) to find the best match (closer to orthogonal) */

        // Vector from the point to its projection on the quad facing the u axis
        const auto point_to_u = projected_on_u - p;

        // Vector from the point to its projection on the quad facing the v axis
        const auto point_to_v = projected_on_v - p;

        // The u-axis of the computed frame
        const auto u = point_to_u.normalized();

        // The v-axis of the computed frame
        auto v = point_to_v.normalized();

        // The w-axis of the computed frame
        const WorldCoordinates w = u.cross(v).normalized();

        // v-axis (recompute the v-axis in case u and v aren't orthogonal
        v = w.cross(u).normalized();

        Matrix<3, 3> m;
        m << u, v, w;

        return m;
    }

    /** Compute the jacobian matrix evaluated at local position {u,v,w}
 * (see interpolation::CanonicalElement::Jacobian for more details).
 * */
    inline
    Matrix<3, 3>
    jacobian (const LocalCoordinates & coordinates) const
    {
        return CanonicalElementType::Jacobian(coordinates, nodes());
    }

    /**
     * Compute the transformation of a local position {u,v,w} to its world position {x,y,z}
     */
    inline
    WorldCoordinates
    T(const LocalCoordinates & coordinates) const
    {
        return CanonicalElementType::interpolate(coordinates, nodes());
    }

    /**
     * Compute an integral approximation by gauss quadrature on the hexahedron of the given evaluation function.
     *
     * @example
     * \code{.cpp}
     * // Integrate the polynomial 1 + 2x + 2xy + 3*z on an hexahedron.
     * float result = Hexahedron(x1, x2, x3, x4, x5, x6, x7, x8).gauss_integrate(
     *   [] (const Hexahedron & hexa, const Hexahedron::LocalCoordinates & coordinates) -> float {
     *     const auto & xi   = coordinates[0];
     *     const auto & eta  = coordinates[1];
     *     const auto & zeta = coordinates[2];
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
     *     ValueType f (const Hexahedron & hexa, const LocalCoordinates & coordinates);
     *
     * Where hexa is a reference to the current hexahadron on which we integrate, and the coordinates u, v and w
     * forms the local position of a sample point on which we want to get the evaluation value.
     *
     * @return The value of the integral computed on this hexahedron.
     *
     */
    template <typename ValueType, typename EvaluateFunctor>
    inline
    ValueType
    gauss_quadrature(EvaluateFunctor f) const
    {
        const auto p0 = MapVector<3>(CanonicalElementType::gauss_nodes[0]);
        const auto w0 = CanonicalElementType::gauss_weights[0];
        const auto detJ0 = jacobian(p0).determinant();
        const auto eval0 = f(*this, p0);
        ValueType result = eval0 * w0 * detJ0;

        for (std::size_t i = 1; i < CanonicalElementType::number_of_gauss_nodes; ++i) {
            const auto p = MapVector<3>(CanonicalElementType::gauss_nodes[i]);
            const auto w = CanonicalElementType::gauss_weights[i];
            const auto detJ = jacobian(p).determinant();
            const auto eval = f(*this, p);
            result += eval * w * detJ;
        }

        return result;
    }

private:
    template <size_t index, typename ...Nodes, REQUIRES(sizeof...(Nodes) >= 1)>
    inline
    void construct_from_nodes(const WorldCoordinates & first_node, Nodes&&...remaining_nodes) {
        p_nodes.row(index) = first_node;
        construct_from_nodes<index+1>(std::forward<Nodes>(remaining_nodes)...);
    }

    template <size_t index>
    inline
    void construct_from_nodes(const WorldCoordinates & last_node) {
        p_nodes.row(index) = last_node;
    }

private:
    Matrix<NumberOfNodes, 3> p_nodes;
};

template <
        typename ...Nodes,
        REQUIRES(8 == sizeof...(Nodes))
>
Hexahedron(Nodes&&...remaining_nodes) -> Hexahedron<interpolation::Hexahedron8>;

Hexahedron() -> Hexahedron<interpolation::Hexahedron8>;

} // namespace caribou::geometry
#endif //CARIBOU_GEOMETRY_HEXAHEDRON_H
