#ifndef CARIBOU_GEOMETRY_INTERNAL_BASEHEXAHEDRON_H
#define CARIBOU_GEOMETRY_INTERNAL_BASEHEXAHEDRON_H

#include <Caribou/config.h>
#include <Caribou/Algebra/Vector.h>

namespace caribou {
namespace geometry {
namespace internal {


template<typename CanonicalElementType, typename HexahedronType>
struct BaseHexahedron : public CanonicalElementType
{
    static constexpr INTEGER_TYPE NumberOfNodes = HexahedronType::NumberOfNodes;
    using CanonicalElement = CanonicalElementType;
    using NodeType = caribou::geometry::Node<3>;
    using QuadType = Quad<3, typename CanonicalElementType::QuadType>;
    using Index = std::size_t ;
    using Real = FLOATING_POINT_TYPE;

    using LocalCoordinates = algebra::Vector<3, Real>;
    using WorldCoordinates = algebra::Vector<3, Real>;


    /**
     * Get the ith quadrangle face.
     */
    inline
    QuadType
    face(Index index) const
    {
        const auto & face_indices = HexahedronType::faces[index];
        std::array<NodeType, QuadType::NumberOfNodes> quad_nodes;
        for (std::size_t i = 0; i < QuadType::NumberOfNodes; ++i)
            quad_nodes[i] = self().node(face_indices[i]);

        return QuadType(quad_nodes);
    }

    /** Compute the volume of the hexahedron */
    inline
    FLOATING_POINT_TYPE
    volume() const
    {
        FLOATING_POINT_TYPE v = 0.;
        for (std::size_t gauss_node_id = 0; gauss_node_id < CanonicalElementType::gauss_nodes.size(); ++gauss_node_id) {
            const auto &gauss_node   = CanonicalElementType::gauss_nodes[gauss_node_id];
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

} // namespace internal
} // namespace geometry
} // namespace caribou
#endif //CARIBOU_GEOMETRY_INTERNAL_BASEHEXAHEDRON_H
