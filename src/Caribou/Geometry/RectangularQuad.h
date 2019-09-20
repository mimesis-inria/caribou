#ifndef CARIBOU_GEOMETRY_RECTANGULARQUAD_H
#define CARIBOU_GEOMETRY_RECTANGULARQUAD_H

#include <Eigen/Core>

#include <Caribou/config.h>
#include <Caribou/Traits.h>
#include <Caribou/Geometry/Interpolation/Quad.h>

namespace caribou::geometry {

template <size_t Dim, typename CanonicalElementType = interpolation::Quad4>
struct RectangularQuad : public CanonicalElementType
{
    static constexpr INTEGER_TYPE NumberOfNodes = CanonicalElementType::NumberOfNodes;

    using LocalCoordinates = typename CanonicalElementType::LocalCoordinates;
    using WorldCoordinates = Eigen::Matrix<FLOATING_POINT_TYPE, Dim, 1>;

    template<int nRows, int nColumns, int Options=Eigen::RowMajor>
    using Matrix = Eigen::Matrix<FLOATING_POINT_TYPE, nRows, nColumns, Options>;

    template<int nRows, int Options=0>
    using Vector = Eigen::Matrix<FLOATING_POINT_TYPE, nRows, 1, Options>;

    static_assert(Dim == 2 or Dim == 3, "Only 2D and 3D quads are supported.");

    using Size = Vector<Dim>;
    using Rot = Matrix <Dim, Dim>;

    constexpr
    RectangularQuad()
        : p_center (WorldCoordinates::Zero()), p_H (Size::Constant(2)), p_R (Rot::Identity())
    {}

    constexpr
    RectangularQuad(const WorldCoordinates & center, const Size & dimensions, const Rot & rotation)
        : p_center (center), p_H (dimensions), p_R (rotation)
    {}

    constexpr
    RectangularQuad(const WorldCoordinates & center, const Size & dimensions)
        : p_center (center), p_H (dimensions), p_R (Rot::Identity())
    {}

    constexpr
    RectangularQuad(const WorldCoordinates & center)
        : p_center (center), p_H (Size::Constant(2)), p_R (Rot::Identity())
    {}

    inline
    WorldCoordinates
    node(UNSIGNED_INTEGER_TYPE index) const
    {
        return T(LocalCoordinates(&CanonicalElementType::nodes[index][0]));
    }

    /** Get a reference to the set of nodes */
    inline
    Matrix<NumberOfNodes, Dim>
    nodes() const
    {
        Matrix<NumberOfNodes, Dim> m;
        for (size_t i = 0; i < CanonicalElementType::NumberOfNodes; ++i)
            m.row(i) = node(i).transpose();
        return m;
    }

    /** Compute the center position **/
    auto
    center() const noexcept
    {
        return T(LocalCoordinates({0., 0.}));
    }

    /**
     * Get the local coordinates frame (a.k.a. the rotation matrix) positioned at the center of the quad
     */
    inline
    Rot
    frame() const
    {
        return p_R;
    }

    /**
     * Get the size (hx, hy, hz) of the quad
     */
    inline
    Size
    H() const
    {
        return p_H;
    }

    /**
     * Compute the transformation of a local position {u} to its world position {x,y,z}
     */
    inline
    WorldCoordinates
    T(const LocalCoordinates & coordinates) const
    {
        return p_center + ((p_R * coordinates).array()*(p_H/2.).array()).matrix();
    }

    /**
     * Compute the inverse transformation of a world position {x,y,z} to its local position {u,v,w}
     */
    inline
    LocalCoordinates
    Tinv(const WorldCoordinates & coordinates) const
    {
        return p_R.transpose() * ((coordinates - p_center).array() / (p_H/2.).array()).matrix();
    }

    /**
     * Returns true if the given world coordinates are within the quad's boundaries, false otherwise.
     */
    inline bool
    contains(const WorldCoordinates & coordinates) const
    {
        const LocalCoordinates c = Tinv(coordinates);
        return IN_CLOSED_INTERVAL(-1, c[0], 1) and
               IN_CLOSED_INTERVAL(-1, c[1], 1);
    }

    /**
     * Compute the jacobian matrix evaluated at local position {u,v}
     *
     * For a rectangular quads, the jacobian matrix is constant and is defined as
     *
     *     1 | hx 0  |
     * J = - | 0  hy |
     *     2
     *
     * where hx, hy are the dimension of the edges 0-1 and 0-3 respectively.
     */
    inline
    Matrix<Dim, 2>
    jacobian (const LocalCoordinates & /*coordinates*/) const
    {
        return jacobian();
    }

    /**
     * Compute the jacobian matrix evaluated at local position {u,v}
     *
     * For a rectangular quads, the jacobian matrix is constant and is defined as
     *
     *     1 | hx 0  |
     * J = - | 0  hy |
     *     2
     *
     * where hx, hy are the dimension of the edges 0-1 and 0-3 respectively.
     */
    inline
    Matrix<Dim, 2>
    jacobian () const {
        return (1 / 2. * p_H).asDiagonal();
    }

private:
    WorldCoordinates p_center; ///< Position of the center point of the quad
    Size p_H; ///< Size of the quad {hx, hy, hz}
    Rot p_R; ///< Rotation matrix (a.k.a. the local coordinates frame) at the center of the quad
};

} // namespace caribou::geometry
#endif //CARIBOU_GEOMETRY_RECTANGULARQUAD_H
