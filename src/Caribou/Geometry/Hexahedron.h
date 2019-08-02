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

    using QuadType = Quad<3, typename CanonicalElementType::QuadType>;

    template<int nRows, int nColumns, int Options=0>
    using Matrix = Eigen::Matrix<FLOATING_POINT_TYPE, nRows, nColumns, Options>;

    template<int nRows, int nColumns>
    using Map = Eigen::Map<const Matrix<nRows, nColumns, Eigen::RowMajor>>;

    template<int nRows, int Options=0>
    using Vector = Eigen::Matrix<FLOATING_POINT_TYPE, nRows, 1, Options>;

    template<int nRows>
    using MapVector = Eigen::Map<const Vector<nRows, Eigen::ColMajor>>;

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

    /** Get the Node at given index */
    inline
    auto
    node(UNSIGNED_INTEGER_TYPE index) const
    {
        return p_nodes.row(index).transpose();
    }

    /** Get the Node at given index */
    inline
    auto
    node(UNSIGNED_INTEGER_TYPE index)
    {
        return p_nodes.row(index).transpose();
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
    QuadType
    face(UNSIGNED_INTEGER_TYPE index) const
    {
        const auto & face_indices = CanonicalElementType::faces[index];

        Matrix<QuadType::NumberOfNodes, 3> m;
        for (std::size_t i = 0; i < QuadType::NumberOfNodes; ++i)
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
    Eigen::Matrix<FLOATING_POINT_TYPE, 3, 3>
    frame() const
    {
        const auto hexa_center = T( LocalCoordinates {0,0,0} ); // Hexahedron's center position

        const auto quad_faced_to_u_axis = face(2); // Quad that lies in front of the u axis
        const auto face_u_center = quad_faced_to_u_axis.T({0,0}); // Center position of this quad

        const auto quad_faced_to_v_axis = face(4); // Quad that lies in front of the v axis
        const auto face_v_center = quad_faced_to_v_axis.T({0,0}); // Center position of this quad

        /* @todo(jnbrunet2000@gmail.com): select between the pairs of axis (center-to-u, center-to-v),
          (center-to-u, center-to-w) and (center-to-v, center-to-w) to find the best match (closer to orthogonal) */

        // Vector from the hexa's center to the center of the quad faced to the u axis
        const auto center_to_u = face_u_center - hexa_center;

        // Vector from the hexa's center to the center of the quad faced to the v axis
        const auto center_to_v = face_v_center - hexa_center;

        // u-axis
        const auto u = center_to_u.normalized();

        // v-axis
        auto v = center_to_v.normalized();

        // w-axis
        const auto w = u.cross(v).normalized();

        // v-axis (recompute the v-axis in case u and v aren't orthogonal
        v = w.cross(u).normalized();

        Eigen::Matrix<FLOATING_POINT_TYPE, 3, 3> m;
        m << u, v, w;

        return m;
    }

    /** Compute the jacobian matrix evaluated at local position {u,v,w}
 * (see interpolation::CanonicalElement::Jacobian for more details).
 * */
    inline
    Eigen::Matrix<FLOATING_POINT_TYPE, 3, 3>
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
        return CanonicalElementType::interpolate_at_local_position(coordinates, nodes());
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
    template <typename EvaluateFunctor>
    inline
    auto
    gauss_quadrature(EvaluateFunctor f) const
    {
        const auto p0 = MapVector<3>(CanonicalElementType::gauss_nodes[0]);
        const auto w0 = CanonicalElementType::gauss_weights[0];
        const auto detJ0 = jacobian(p0).determinant();
        const auto eval0 = f(*this, p0);
        auto result = eval0 * w0 * detJ0;

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
