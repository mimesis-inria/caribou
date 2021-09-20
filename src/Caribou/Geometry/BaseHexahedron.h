#pragma once

#include <Caribou/config.h>
#include <Caribou/constants.h>
#include <Caribou/Geometry/Element.h>
#include <Caribou/Geometry/Triangle.h>
#include <Eigen/Core>
#include <Eigen/Dense>

namespace caribou::geometry {

/**
 * Base class for hexahedral elements.
 * @tparam Derived Derived hexahedral type (usually a leaf type, which is a type that is used by the user, eg Hexahedron).
 */
template<typename Derived>
struct BaseHexahedron : public Element<Derived> {
    // Types
    using Base = Element<Derived>;

    using LocalCoordinates = typename Base::LocalCoordinates;
    using WorldCoordinates = typename Base::WorldCoordinates;

    using GaussNode = typename Base::GaussNode;
    using Mat3x3   = Eigen::Matrix<double, 3, 3>;

    template <UNSIGNED_INTEGER_TYPE Dim>
    using Vector = typename Base::template Vector<Dim>;

    template <UNSIGNED_INTEGER_TYPE Rows, UNSIGNED_INTEGER_TYPE Cols>
    using Matrix = typename Base::template Matrix<Rows, Cols>;

    // Constants
    static constexpr auto CanonicalDimension = Base::CanonicalDimension;
    static constexpr auto Dimension = Base::Dimension;
    static constexpr auto NumberOfNodesAtCompileTime = Base::NumberOfNodesAtCompileTime;
    static constexpr auto NumberOfGaussNodesAtCompileTime = Base::NumberOfGaussNodesAtCompileTime;

    static_assert(Dimension == 3, "Hexahedrons can only be of dimension 3.");

    /** Default empty constructor */
    BaseHexahedron() = default;

    /** Constructor from an Eigen matrix containing the positions of the hexahedron's nodes */
    template<typename EigenType, REQUIRES(EigenType::RowsAtCompileTime == NumberOfNodesAtCompileTime)>
    explicit BaseHexahedron(Eigen::EigenBase<EigenType> & nodes) :p_nodes(nodes.derived().template cast<typename Base::Scalar>()) {}

    /** Constructor from an Eigen matrix containing the positions of the hexahedron's nodes */
    template<typename EigenType, REQUIRES(EigenType::RowsAtCompileTime == NumberOfNodesAtCompileTime)>
    explicit BaseHexahedron(const Eigen::EigenBase<EigenType> & nodes) :p_nodes(nodes.derived().template cast<typename Base::Scalar>()) {}

    /** Constructor from a serie of nodes. */
    template <
        typename ...Nodes,
        REQUIRES(NumberOfNodesAtCompileTime == sizeof...(Nodes)+1)
    >
    explicit BaseHexahedron(const WorldCoordinates & first_node, Nodes&&...remaining_nodes)
    {
        construct_from_nodes<0>(first_node, std::forward<Nodes>(remaining_nodes)...);
    }

    // Public methods common to all hexahedral types

    /**
     * Get the list of node indices of the faces.
     * \sa Element::boundary_elements_node_indices
     */
    inline auto faces() const {
        return self().get_boundary_elements_nodes();
    }

    /**
     * Extract the orthogonal frame of the element by computing the cross product of the unit vectors
     * from the center position to its projection on opposite faces.
     *
     * \sa BaseHexahedron::frame(const LocalCoordinates & local_point)
     *
     * \warning If the hexahedron isn't rectangular, the frame extracted by this function will be a rough
     *          approximation that could be far from the real solution, especially for strongly deformed or
     *          inverted hexahedrons.
     */
    inline auto frame() const -> Matrix<3, 3>
    {
        return frame(LocalCoordinates::Zeros());
    }

    /**
     * Extract the frame positioned at the given position (in local coordinates) on the hexa by computing the cross
     * product of the unit vectors from the given position its projection on opposite faces.
     *
     * This function will return a matrix of the form:
     * | ux vx wx |
     * | uy vy wy |
     * | uz vz wz |
     *
     * Where (ux, uy, uz), (vx, vy, vz) and (wx, wy, wz) are orthogonal unitary vectors representing
     * the u, v and w frames in the current hexa. If the hexa is rectangular and not rotated, this matrix is the
     * Identity matrix. If it is rectangular but rotated, rotating the hexa by the transposed of this frame should
     * align the u,v ,w axis to the x,y,z world frame (identity matrix).
     *
     * \warning If the hexahedron isn't rectangular, the frame extracted by this function will be a rough
     *          approximation that could be far from the real solution, especially for strongly deformed or
     *          inverted hexahedrons.
     */
    inline auto frame(const LocalCoordinates & local_point) const -> Matrix<3, 3>
    {
        // Position of the point inside the hexahedron where the frame should be computed
        const auto p = this->world_coordinates( local_point );

        // Project of the point on the quad facing the u axis
        const auto projected_on_u = this->world_coordinates({1, local_point[1],local_point[2]});

        // Project of the point on the quad facing the v axis
        const auto projected_on_v = this->world_coordinates({local_point[0], 1,local_point[2]});

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

        /* Eigen::HouseholderQR<Matrix<3, 3>> qr(m.rows(), m.cols());
        qr.compute(m);

        Matrix<3, 3> Q = qr.householderQ(); */
        
        /* Matrix<3, 3> R = qr.householderR();
        Matrix<3, 3> t_Q = Q.transpose();
        Matrix<3, 3> res = Q*t_Q;
        std::cout << "This is the original Matrix: \n" << m << std::endl;
        std::cout << "This is its Orthogonal component: \n" << Q << std::endl;
        std::cout << "The dot product: \n" << res << std::endl; */

        return m;
    }






    /* inline auto get_local_base() -> Mat3x3 {
        const auto ex = this->world_coordinates({1, 0, 0}); 
        const auto ey = this->world_coordinates({0, 1, 0}); 
        const auto ez = this->world_coordinates({0, 0, 1});

        Mat3x3 base; 
        base << ex, ey, ez;
        return base; 
    } */

    /**
     * Test if the cube intersects the given 3D segment (in world coordinates)
     *
     * @note  based on polygon_intersects_cube by Don Hatch (January 1994)
     */
    [[nodiscard]]
    inline auto intersects(const Segment<_3D, Linear> & segment, const FLOATING_POINT_TYPE eps=0) const -> bool {
        return intersects_local_segment(self().local_coordinates(segment.node(0)), self().local_coordinates(segment.node(1)), eps);
    }

    /**
     * Test if the cube intersects the given 3D triangle (in world coordinates)
     *
     * @note  based on polygon_intersects_cube by Don Hatch (January 1994)
     */
    [[nodiscard]]
    inline auto intersects(const Triangle<_3D, Linear> & t, const FLOATING_POINT_TYPE eps=0) const -> bool {
        LocalCoordinates nodes[3];
        for (UNSIGNED_INTEGER_TYPE i = 0; i < 3; ++i) {
            nodes[i] = self().local_coordinates(t.node(i));
        }

        auto normal = (nodes[1] - nodes[2]).cross(nodes[2] - nodes[0]).eval().normalized();

        return intersects_local_polygon<3>(nodes, normal, eps);
    }

private:
    // Implementations
    friend struct Element<Derived>;
    [[nodiscard]]
    inline auto get_number_of_nodes() const {return NumberOfNodesAtCompileTime;}
    inline auto get_number_of_gauss_nodes() const {return NumberOfGaussNodesAtCompileTime;}
    inline auto get_node(const UNSIGNED_INTEGER_TYPE & index) const {return WorldCoordinates(p_nodes.row(index));};
    inline auto get_nodes() const -> const auto & {return p_nodes;};
    inline auto get_center() const {return Base::world_coordinates(LocalCoordinates(0, 0, 0));};
    inline auto get_number_of_boundary_elements() const -> UNSIGNED_INTEGER_TYPE {return 6;};
    inline auto get_contains_local(const LocalCoordinates & xi, const FLOATING_POINT_TYPE & eps) const -> bool {
        const auto & u = xi[0];
        const auto & v = xi[1];
        const auto & w = xi[2];
        return IN_CLOSED_INTERVAL(-1-eps, u, 1+eps) and
               IN_CLOSED_INTERVAL(-1-eps, v, 1+eps) and
               IN_CLOSED_INTERVAL(-1-eps, w, 1+eps);
    }

    /**
    * Test if the cube intersects the given 3D segment (in the hexahedron's local coordinates).
    *
    * @note  based on polygon_intersects_cube by Don Hatch (January 1994)
    */
    [[nodiscard]]
    inline auto intersects_local_segment(LocalCoordinates v0, LocalCoordinates v1, const FLOATING_POINT_TYPE & eps) const -> bool;

    /**
     * Test if the cube intersects the given 3D polygon (in world coordinates).
     *
     * @note  based on polygon_intersects_cube by Don Hatch (January 1994)
     * @tparam NNodes Number of nodes in the polygon
     * @param nodes The nodes of the polygon
     * @param polynormal Vector perpendicular to the polygon.  It need not be of unit length.
     * @return True if the polygon intersects the cube, false otherwise.
     */
    template <int NNodes>
    inline auto intersects_local_polygon(const LocalCoordinates nodes[NNodes], const Vector<3> & polynormal, const FLOATING_POINT_TYPE & eps) const -> bool;

    auto self() -> Derived& { return *static_cast<Derived*>(this); }
    auto self() const -> const Derived& { return *static_cast<const Derived*>(this); }

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
protected:
    Matrix<NumberOfNodesAtCompileTime, Dimension> p_nodes;
};

template<typename Derived>
auto BaseHexahedron<Derived>::intersects_local_segment(LocalCoordinates v0, LocalCoordinates v1, const FLOATING_POINT_TYPE & eps) const -> bool
{
    const auto edge = (v1 - v0);
    INTEGER_TYPE edge_signs[3];

    for (UNSIGNED_INTEGER_TYPE i = 0; i < 3; ++i) {
        edge_signs[i] = (edge[i] < 0) ? -1 : 1;
    }

    for (UNSIGNED_INTEGER_TYPE i = 0; i < 3; ++i) {

        if (v0[i] * edge_signs[i] >  1+eps) return false;
        if (v1[i] * edge_signs[i] < -1-eps) return false;
    }


    for (UNSIGNED_INTEGER_TYPE i = 0; i < 3; ++i) {
        FLOATING_POINT_TYPE rhomb_normal_dot_v0, rhomb_normal_dot_cubedge;

        const UNSIGNED_INTEGER_TYPE iplus1 = (i + 1) % 3;
        const UNSIGNED_INTEGER_TYPE iplus2 = (i + 2) % 3;

        rhomb_normal_dot_v0 = edge[iplus2] * v0[iplus1] -
                              edge[iplus1] * v0[iplus2];

        rhomb_normal_dot_cubedge = edge[iplus2] * edge_signs[iplus1] +
                                   edge[iplus1] * edge_signs[iplus2];

        const auto r = (rhomb_normal_dot_v0*rhomb_normal_dot_v0) - (rhomb_normal_dot_cubedge*rhomb_normal_dot_cubedge);
        if (r > eps)
            return false;
    }

    return true;
}

#define seg_contains_point(a,b,x) (((b)>(x)) - ((a)>(x)))

template<typename Derived>
template <int NNodes>
inline auto BaseHexahedron<Derived>::intersects_local_polygon(const LocalCoordinates nodes[NNodes], const Vector<3> & polynormal, const FLOATING_POINT_TYPE & eps) const -> bool {

    // Check if any edges of the polygon intersect the hexa
    for (UNSIGNED_INTEGER_TYPE i = 0; i < NNodes; ++i) {
        const auto & p1 = nodes[i];
        const auto & p2 = nodes[(i+1)%NNodes];
        if (intersects_local_segment(p1, p2, eps))
            return true;
    }

    // Check that if the polygon's plane intersect the cube diagonal that is the closest of being perpendicular to the
    // plane of the polygon.
    Vector<3> best_diagonal;
    best_diagonal << (((polynormal[0]) < 0) ? -1 : 1),
            (((polynormal[1]) < 0) ? -1 : 1),
            (((polynormal[2]) < 0) ? -1 : 1);

    // Check if the intersection point between the two planes lies inside the cube
    const FLOATING_POINT_TYPE t = polynormal.dot(nodes[0]) / polynormal.dot(best_diagonal);
    if (!IN_CLOSED_INTERVAL(-1-eps, t, 1+eps))
        return false;

    // Check if the intersection point between the two planes lies inside the polygon
    LocalCoordinates p = best_diagonal * t;
    const LocalCoordinates abspolynormal = polynormal.array().abs();
    int zaxis, xaxis, yaxis;
    if (abspolynormal[0] > abspolynormal[1])
        zaxis = (abspolynormal[0] > abspolynormal[2]) ? 0 : 2;
    else
        zaxis = (abspolynormal[1] > abspolynormal[2]) ? 1 : 2;

    if (polynormal[zaxis] < 0) {
        xaxis = (zaxis+2)%3;
        yaxis = (zaxis+1)%3;
    }
    else {
        xaxis = (zaxis+1)%3;
        yaxis = (zaxis+2)%3;
    }

    int count = 0;
    for (UNSIGNED_INTEGER_TYPE i = 0; i < NNodes; ++i) {
        const auto & p1 = nodes[i];
        const auto & p2 = nodes[(i+1)%NNodes];

        if (const int xdirection = seg_contains_point(p1[xaxis]-eps, p2[xaxis]+eps, p[xaxis]))
        {
            if (seg_contains_point(p1[yaxis]-eps, p2[yaxis]+eps, p[yaxis]))
            {
                if (xdirection * (p[xaxis]-p1[xaxis])*(p2[yaxis]-p1[yaxis]) <=
                    xdirection * (p[yaxis]-p1[yaxis])*(p2[xaxis]-p1[xaxis]))
                    count += xdirection;
            }
            else
            {
                if (p2[yaxis] <= p[yaxis])
                    count += xdirection;
            }
        }

    }

    return (count != 0);
}

#undef seg_contains_point

}