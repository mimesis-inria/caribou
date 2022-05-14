#pragma once

#include <Caribou/config.h>
#include <Caribou/constants.h>
#include <Caribou/Geometry/Element.h>
#include <Caribou/Geometry/Triangle.h>
#include <Eigen/Core>


namespace caribou::geometry {

/**
 * Base class for rectangular hexahedral elements.
 * @tparam Derived Derived rectangular hexahedral type (usually a leaf type, which is a type that is used by the user,
 * eg RectangularHexahedron).
 */
template<typename Derived>
struct BaseRectangularHexahedron : public Element<Derived> {
    // Types
    using Base = Element<Derived>;

    using LocalCoordinates = typename Base::LocalCoordinates;
    using WorldCoordinates = typename Base::WorldCoordinates;

    using GaussNode = typename Base::GaussNode;

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

    using Size = Vector<3>;
    using Rotation = Matrix<3,3>;

    /** Default empty constructor */
    BaseRectangularHexahedron() : p_center(0, 0, 0), p_H(2,2,2), p_R(Rotation::Identity()) {}

    /** Constructor by specifying the center point */
    explicit BaseRectangularHexahedron(WorldCoordinates center) : p_center(center), p_H(Size::Constant(2)), p_R(Rotation::Identity()) {}

    /** Constructor by specifying the center point and the size (hx, hy, hz) */
    BaseRectangularHexahedron(WorldCoordinates center, Size H) : p_center(center), p_H(H), p_R(Rotation::Identity()) {}

    /** Constructor by specifying the center point and the rotation */
    BaseRectangularHexahedron(WorldCoordinates center, Rotation R) : p_center(center), p_H(Size::Constant(2)), p_R(R) {}

    /** Constructor by specifying the center point, the size (hx, hy, hz) and the rotation */
    BaseRectangularHexahedron(WorldCoordinates center, Size H, Rotation R) : p_center(center), p_H(H), p_R(R) {}

    // Public methods common to all rectangular hexahedral types

    /**
     * Get the list of node indices of the faces.
     * \sa Element::boundary_elements_node_indices
     */
    inline auto faces() const {
        return self().boundary_elements_node_indices();
    }

    /** Get the rotation frame of the quad */
    inline auto rotation() const -> const Rotation  & {
        return p_R;
    };

    /** Get the size (hx, hy) of the quad */
    inline auto size() const -> const Size  & {
        return p_H;
    };

    /** Get the world coordinates of a point from its local coordinates. */
    inline auto world_coordinates(const LocalCoordinates & coordinates) const -> WorldCoordinates {
        return p_center + p_R * (coordinates.cwiseProduct(p_H / 2.));
    }

    /** Get the local coordinates of a point from its world coordinates. */
    inline auto local_coordinates(const WorldCoordinates & coordinates) const -> LocalCoordinates {
        return (p_R.transpose()*(coordinates - p_center)).cwiseQuotient(p_H / 2.);
    }

    /** Returns true if the given world coordinates are within the hexahedron's boundaries, false otherwise. */
    inline auto contains(const WorldCoordinates & coordinates) const -> bool
    {
        const auto c = local_coordinates(coordinates);
        return self().contains_local(c);
    }

    /**
     * Test if the cube intersects the given 3D segment (in world coordinates)
     *
     * @note  based on polygon_intersects_cube by Don Hatch (January 1994)
     */
    [[nodiscard]]
    inline auto intersects(const Segment<_3D> & segment, const FLOATING_POINT_TYPE eps=0) const -> bool {
        return intersects_local_segment(local_coordinates(segment.node(0)), local_coordinates(segment.node(1)), eps);
    }

    /**
     * Test if the cube intersects the given 3D triangle (in world coordinates)
     *
     * @note  based on polygon_intersects_cube by Don Hatch (January 1994)
     */
    [[nodiscard]]
    inline auto intersects(const Triangle<_3D> & t, const FLOATING_POINT_TYPE eps=0) const -> bool {
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
    [[nodiscard]]
    inline auto get_number_of_gauss_nodes() const {return NumberOfGaussNodesAtCompileTime;}
    inline auto get_center() const {return p_center;};
    [[nodiscard]]
    inline auto get_number_of_boundary_elements() const -> UNSIGNED_INTEGER_TYPE {return 6;};
    inline auto get_contains_local(const LocalCoordinates & xi, const FLOATING_POINT_TYPE & eps) const -> bool {
        const auto & u = xi[0];
        const auto & v = xi[1];
        const auto & w = xi[2];
        return IN_CLOSED_INTERVAL(-1-eps, u, 1+eps) and
               IN_CLOSED_INTERVAL(-1-eps, v, 1+eps) and
               IN_CLOSED_INTERVAL(-1-eps, w, 1+eps);
    }

    auto self() -> Derived& { return *static_cast<Derived*>(this); }
    auto self() const -> const Derived& { return *static_cast<const Derived*>(this); }

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

protected:
    WorldCoordinates p_center; ///< Position of the center point of the hexahedron
    Size p_H; ///< Size of the hexahedron {hx, hy, hz}
    Rotation p_R; ///< Rotation matrix (a.k.a. the local coordinates frame) at the center of the hexahedron
};

template<typename Derived>
auto BaseRectangularHexahedron<Derived>::intersects_local_segment(LocalCoordinates v0, LocalCoordinates v1, const FLOATING_POINT_TYPE & eps) const -> bool
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
inline auto BaseRectangularHexahedron<Derived>::intersects_local_polygon(const LocalCoordinates nodes[NNodes], const Vector<3> & polynormal, const FLOATING_POINT_TYPE & eps) const -> bool {

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