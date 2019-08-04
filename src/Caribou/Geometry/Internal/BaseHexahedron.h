#ifndef CARIBOU_GEOMETRY_INTERNAL_BASEHEXAHEDRON_H
#define CARIBOU_GEOMETRY_INTERNAL_BASEHEXAHEDRON_H

#include <Eigen/Core>

#include <Caribou/config.h>

namespace caribou::geometry::internal {

template<typename CanonicalElementType, typename HexahedronType>
struct BaseHexahedron : public CanonicalElementType
{
    static constexpr INTEGER_TYPE NumberOfNodes = HexahedronType::NumberOfNodes;
    using CanonicalElement = CanonicalElementType;

    using LocalCoordinates = typename CanonicalElement::LocalCoordinates;
    using WorldCoordinates = Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1>;

    template<int nRows, int Options=0>
    using Vector = Eigen::Matrix<FLOATING_POINT_TYPE, nRows, 1, Options>;

    template<int nRows>
    using MapVector = Eigen::Map<const Vector<nRows, Eigen::ColMajor>>;

    /** Compute the volume of the hexahedron */
    inline
    FLOATING_POINT_TYPE
    volume() const
    {
        FLOATING_POINT_TYPE v = 0.;
        for (std::size_t gauss_node_id = 0; gauss_node_id < CanonicalElementType::number_of_gauss_nodes; ++gauss_node_id) {
            const auto &gauss_node = MapVector<3>(CanonicalElementType::gauss_nodes[gauss_node_id]);
            const auto &gauss_weight = CanonicalElementType::gauss_weights[gauss_node_id];
            const auto J = self().jacobian(gauss_node);
            const auto detJ = J.determinant();

            v += detJ * gauss_weight;
        }
        return v;
    }

private:
    const HexahedronType &self () const
    {
        return static_cast<const HexahedronType &>(*this);
    }
};

} // namespace caribou::geometry::internal
#endif //CARIBOU_GEOMETRY_INTERNAL_BASEHEXAHEDRON_H
