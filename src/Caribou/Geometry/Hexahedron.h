#ifndef CARIBOU_GEOMETRY_HEXAHEDRON_H
#define CARIBOU_GEOMETRY_HEXAHEDRON_H

#include <array>

namespace caribou
{
namespace geometry
{

/**
 * The base hexahedron class (of n nodes) in 3D space.
 *
 * This class will be overloaded with the number of nodes (8-nodes linear, 20-nodes quadratic, etc.) and the
 * type (regular or non-regular).
 */
template <size_t NNodes, typename HexahedronTraits_>
struct Hexahedron
{

    static constexpr size_t NumberOfNodes = NNodes;

    using HexahedronTraits = HexahedronTraits_;
    using PointType = typename HexahedronTraits::PointType;
    using VectorType = typename PointType::VectorType;
    using ValueType = typename VectorType::ValueType;
    using SegmentType = typename HexahedronTraits::SegmentType;

    static_assert(NumberOfNodes >= 8, "A hexahedron must have at least eight nodes.");

    /** Delete the default constructor as it should be overloaded by base class on the number of nodes **/
    Hexahedron () = delete;

    /**
     * Constructor by a list initializer of PointType.
     */
    Hexahedron (const std::initializer_list<PointType> & il) : HexahedronTraits(il)
    {
    }

    /** Constructor by an array of points **/
    Hexahedron (const std::array<PointType, NumberOfNodes> & n) : HexahedronTraits(n)
    {
    }

//    Hexahedron (const std::initializer_list<VectorType> & il)
//    {
//        auto v = std::begin(il);
//        size_t i = 0;
//        for (; v != std::end(il); ++i, ++v) {
//            nodes[i][0] = (*v)[0];
//            nodes[i][1] = (*v)[1];
//            nodes[i][2] = (*v)[2];
//        }
//    }

    /** Scale the hexahedron by s (sx, sy, sz) from the origin **/
    HexahedronTraits
    inline scale(VectorType s) const
    {
        return HexahedronTraits::scale(s);
    }
};

} // namespace geometry

} // namespace caribou

#endif //CARIBOU_GEOMETRY_HEXAHEDRON_H
