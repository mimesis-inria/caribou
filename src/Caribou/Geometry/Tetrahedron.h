#ifndef CARIBOU_GEOMETRY_TETRAHEDRON_H
#define CARIBOU_GEOMETRY_TETRAHEDRON_H

#include <Eigen/Dense>
#include <Caribou/macros.h>
#include <Caribou/config.h>
#include <Caribou/Geometry/Triangle.h>
#include <Caribou/Geometry/Interpolation/Tetrahedron.h>
#include <Caribou/Geometry/Internal/BaseTetrahedron.h>

namespace caribou::geometry {

template <typename CanonicalElementType>
struct Tetrahedron : public internal::BaseTetrahedron<CanonicalElementType::NumberOfNodes, CanonicalElementType, Tetrahedron<CanonicalElementType>>
{
};

template <>
struct Tetrahedron<interpolation::Tetrahedron4> : public internal::BaseTetrahedron<4, interpolation::Tetrahedron4, Tetrahedron<interpolation::Tetrahedron4>>
{
    static constexpr INTEGER_TYPE NumberOfNodes = 4;
    using CanonicalElementType = interpolation::Tetrahedron4;
    using Base = internal::BaseTetrahedron<NumberOfNodes, CanonicalElementType, Tetrahedron<CanonicalElementType>>;
    using WorldCoordinates = typename Base::WorldCoordinates;
    using FaceType = Triangle<3, typename CanonicalElementType::FaceType>;

    using Base::Base;

    /**
     * Compute the inverse transformation of a world position {x,y,z} to its local position {u,v,w}
     */
    inline
    LocalCoordinates
    Tinv(const WorldCoordinates & coordinates) const
    {
        const auto n0 = node(0);
        const auto n1 = node(1);
        const auto n2 = node(2);
        const auto n3 = node(3);

        Matrix<3,3> At;
        At.row(0) = n1-n0;
        At.row(1) = n2-n0;
        At.row(2) = n3-n0;

        Matrix<3,3> A = At.transpose();
        Vector<3>   B = coordinates - n0;

        return A.inverse()*B;
    }
};

template <>
struct Tetrahedron<interpolation::Tetrahedron10> : public internal::BaseTetrahedron<10, interpolation::Tetrahedron10, Tetrahedron<interpolation::Tetrahedron10>>
{
    static constexpr INTEGER_TYPE NumberOfNodes = 10;
    using CanonicalElementType = interpolation::Tetrahedron10;
    using Base = internal::BaseTetrahedron<NumberOfNodes, interpolation::Tetrahedron10, Tetrahedron<CanonicalElementType>>;
    using LocalCoordinates = typename Base::LocalCoordinates;
    using WorldCoordinates = typename Base::WorldCoordinates;
    using FaceType = Triangle<3, typename CanonicalElementType::FaceType>;

    using Base::Base;

    Tetrahedron(const Matrix<10, 3> & m) : Base(m) {}

    template <typename Derived, REQUIRES(Derived::RowsAtCompileTime == 4), REQUIRES(Derived::ColsAtCompileTime == 3)>
    Tetrahedron(const Eigen::MatrixBase<Derived> & m) {
        for (UNSIGNED_INTEGER_TYPE i = 0; i < 4; ++i) {
            p_nodes.row(i) = m.row(i);
        }

        Tetrahedron<interpolation::Tetrahedron4> linear_tetra(m);

        for (UNSIGNED_INTEGER_TYPE i = 4; i < 10; ++i) {
            p_nodes.row(i) = linear_tetra.T(LocalCoordinates(&(interpolation::Tetrahedron10::nodes[i][0])));
        }
    }

    inline bool can_be_converted_to_linear() const {
        for (const auto & edge : CanonicalElementType::edges) {
            const WorldCoordinates edge0 = (Base::node(edge[2]) - Base::node(edge[0])).normalized();
            const WorldCoordinates edge1 = (Base::node(edge[1]) - Base::node(edge[0])).normalized();

            const FLOATING_POINT_TYPE d = (edge0[0]*edge1[0] + edge0[1]*edge1[1] + edge0[2]*edge1[2]) - 1.;
            if (not IN_CLOSED_INTERVAL(-1e-10, d, 1e-10)) {
                return false;
            }
        }

        return true;
    }
};

} // namespace caribou::geometry

#endif //CARIBOU_GEOMETRY_TETRAHEDRON_H
