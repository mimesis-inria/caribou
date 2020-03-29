#pragma once

#include <Caribou/config.h>
#include <Caribou/Constants.h>
#include <Caribou/Geometry/Element.h>
#include <Eigen/Core>


namespace caribou::geometry {

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
        return p_center + p_R * (coordinates.cwiseProduct(p_H / 2.));
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
protected:
    WorldCoordinates p_center; ///< Position of the center point of the hexahedron
    Size p_H; ///< Size of the hexahedron {hx, hy, hz}
    Rotation p_R; ///< Rotation matrix (a.k.a. the local coordinates frame) at the center of the hexahedron
};

}