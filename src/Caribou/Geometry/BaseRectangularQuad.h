#pragma once

#include <Caribou/config.h>
#include <Caribou/Constants.h>
#include <Caribou/Geometry/Element.h>
#include <Eigen/Core>


namespace caribou::geometry {

template<typename Derived>
struct BaseRectangularQuad : public Element<Derived> {
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

    static_assert(Dimension == 2 or Dimension == 3, "Quads can only be of dimension 2 or 3.");

    using Size = Vector<CanonicalDimension>;
    using Rotation = Matrix<Dimension,Dimension>;

    /** Default empty constructor */
    BaseRectangularQuad() : p_center(WorldCoordinates::Constant(0)), p_H(Size::Constant(2)), p_R(Rotation::Identity()) {}

    /** Constructor by specifying the center point */
    explicit BaseRectangularQuad(WorldCoordinates center) : p_center(center), p_H(Size::Constant(2)), p_R(Rotation::Identity()) {}

    /** Constructor by specifying the center point and the size (hx, hy, hz) */
    BaseRectangularQuad(WorldCoordinates center, Size H) : p_center(center), p_H(H), p_R(Rotation::Identity()) {}

    /** Constructor by specifying the center point and the rotation */
    BaseRectangularQuad(WorldCoordinates center, Rotation R) : p_center(center), p_H(Size::Constant(2)), p_R(R) {}

    /** Constructor by specifying the center point, the size (hx, hy, hz) and the rotation */
    BaseRectangularQuad(WorldCoordinates center, Size H, Rotation R) : p_center(center), p_H(H), p_R(R) {}

    // Public methods

    /** Get the rotation frame of the quad */
    inline auto rotation() const -> const Rotation  & {
        return p_R;
    };

    /** Get the size (hx, hy) of the quad */
    inline auto size() const -> const Size  & {
        return p_H;
    };

    /** Compute the transformation of a local position {u,v,w} to its world position {x,y,z} */
    inline auto T(const LocalCoordinates & coordinates) const -> WorldCoordinates {
        if constexpr (Dimension == 2) {
            return p_center + p_R * (coordinates.cwiseProduct(p_H / 2.));
        } else {
            WorldCoordinates p;
            p.template block<CanonicalDimension, 1>(0,0) = coordinates.cwiseProduct(p_H / 2.);
            p[2] = 0.;
            return p_center + p_R * p;
        }
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
    inline auto get_number_of_boundary_elements() const -> UNSIGNED_INTEGER_TYPE {return 4;};
protected:
    WorldCoordinates p_center; ///< Position of the center point of the quad
    Size p_H; ///< Size of the quad {hx, hy}
    Rotation p_R; ///< Rotation matrix (a.k.a. the local coordinates frame) at the center of the quad
};
}