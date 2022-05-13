#pragma once

#include <Caribou/config.h>
#include <Caribou/Geometry/Base/BaseSegment.h>
#include <Eigen/Core>

namespace caribou::geometry {

template<UNSIGNED_INTEGER_TYPE _Dimension>
struct Segment;

template<UNSIGNED_INTEGER_TYPE _Dimension>
 struct traits<Segment <_Dimension>> {
     static constexpr UNSIGNED_INTEGER_TYPE CanonicalDimension = 1;
     static constexpr UNSIGNED_INTEGER_TYPE Dimension = _Dimension;
     static constexpr INTEGER_TYPE NumberOfNodesAtCompileTime = 2;
     static constexpr INTEGER_TYPE NumberOfGaussNodesAtCompileTime = 1;
 };

/**
 * Linear segment
 *
 * \verbatim
 * P1 : 0-----+-----1 --> u
 * \endverbatim
 *
 * @tparam _Dimension The world coordinates dimension
 */
template<UNSIGNED_INTEGER_TYPE _Dimension>
struct Segment: public BaseSegment<Segment <_Dimension>> {
    // Types
    using Base = BaseSegment<Segment <_Dimension>>;
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

    // Constructors
    using Base::Base;
    Segment() : Base() {
        if constexpr (Dimension == 1) {
            this->p_nodes[0] = -1;
            this->p_nodes[1] = +1;
        } else if constexpr (Dimension == 2) {
            this->p_nodes.row(0) = WorldCoordinates(-1, 0);
            this->p_nodes.row(1) = WorldCoordinates(+1, 0);
        } else {
            this->p_nodes.row(0) = WorldCoordinates(-1, 0, 0);
            this->p_nodes.row(1) = WorldCoordinates(+1, 0, 0);
        }
    };


private:
    // Implementations
    friend struct Element<Segment <_Dimension>>;
    friend struct BaseSegment<Segment <_Dimension>>;
    inline auto get_L(const LocalCoordinates & xi) const -> Vector<NumberOfNodesAtCompileTime> {
        const auto  & u = xi[0];
        return {
            static_cast<FLOATING_POINT_TYPE>(1/2. * (1 - u)),
            static_cast<FLOATING_POINT_TYPE>(1/2. * (1 + u))
        };
    };

    inline auto get_dL(const LocalCoordinates & /*xi*/) const -> Matrix<NumberOfNodesAtCompileTime, CanonicalDimension> {
        return {
            static_cast<FLOATING_POINT_TYPE>(-1/2.),
            static_cast<FLOATING_POINT_TYPE>(+1/2.)
        };
    };

    inline auto get_gauss_nodes() const -> const auto & {
        static std::vector<GaussNode> gauss_nodes {
            GaussNode {LocalCoordinates(0), 2} // Node 0
        };
        return gauss_nodes;
    }
};

}
