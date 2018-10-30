#ifndef CARIBOU_GEOMETRY_HEXAHEDRON_H
#define CARIBOU_GEOMETRY_HEXAHEDRON_H

#include <Caribou/Geometry/Point.h>

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
template <size_t NNodes, typename TVector>
struct Hexahedron
{

    static constexpr size_t NumberOfNodes = NNodes;

    using VectorType = TVector;
    using PointType = Point<3, VectorType>;

    static_assert(NumberOfNodes >= 8, "A hexahedron must have at least eight nodes.");

    /** Delete the default constructor as it should be overloaded by base class on the number of nodes **/
    Hexahedron () = delete;

    /**
     * Constructor by a list initializer of PointType.
     *
     * Example:
     *
     * \code{.cpp}
     * Point3D p1, p2, p3, p4, p5, p6, p7, p8;
     * BaseHexahedron<8, false, Point3D::VectorType> hexahedron ({p1, p2, p3, p4, p5, p6, p7, p8});
     * \endcode
     *
     */
    Hexahedron (const std::initializer_list<PointType> & il)
    {
        std::copy(std::begin(il), std::end(il), std::begin(nodes));
    }

    Hexahedron (const std::initializer_list<VectorType> & il)
    {
        auto v = std::begin(il);
        size_t i = 0;
        for (; v != std::end(il); ++i, ++v) {
            nodes[i][0] = (*v)[0];
            nodes[i][1] = (*v)[1];
            nodes[i][2] = (*v)[2];
        }
    }

protected:
    std::array<PointType, NumberOfNodes> nodes;
};

/* Linear hexahedron
 *    7-------6
 *   /|      /|
 *  / |     / |
 * 3--|----2  |
 * |  4----|--5
 * | /     | /
 * 0-------1
 */


template <typename TVector>
struct LinearHexahedron : public Hexahedron<8, TVector>
{
    using Base = Hexahedron<8, TVector>;
//    using Hexahedron<8, TVector>::Hexahedron; // Import the base constructors
    using VectorType = TVector;
    using PointType = Point<3, VectorType>;

    LinearHexahedron (const std::initializer_list<PointType> & il) : Base(il)
    {
    }

    LinearHexahedron (const std::initializer_list<VectorType> & il) : Base(il)
    {
    }

    /** Get the node at index **/
    PointType node(size_t index) const {
        return Base::nodes[index];
    }
};

/**
 * Regular linear hexahedron
 * @tparam TVector
 */
template <typename TVector>
struct RegularLinearHexahedron : public LinearHexahedron<TVector>
{
    using LinearHexahedron<TVector>::LinearHexahedron; // Import the base constructors
};

} // namespace geometry

} // namespace caribou

#endif //CARIBOU_GEOMETRY_HEXAHEDRON_H
