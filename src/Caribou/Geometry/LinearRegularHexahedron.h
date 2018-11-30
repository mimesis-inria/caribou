#ifndef CARIBOU_GEOMETRY_LINEARREGULARHEXAHEDRON_H
#define CARIBOU_GEOMETRY_LINEARREGULARHEXAHEDRON_H

#include <Caribou/Algebra/Matrix.h>
#include <Caribou/Geometry/Point.h>
#include <Caribou/Geometry/Segment.h>
#include <Caribou/Geometry/Quad.h>

#include <Caribou/Geometry/LinearHexahedron.h>

namespace caribou
{
namespace geometry
{

/* Regular Linear hexahedron
 *    7-------6
 *   /|      /|
 *  / |     / |
 * 3--|----2  |
 * |  4----|--5
 * | /     | /
 * 0-------1
 */

struct LinearRegularHexahedron : public LinearHexahedron
{
    using PointType = Point3D;
    using SegmentType = Segment<3>;
    using FaceType = Quad<3>;

    using VectorType= typename PointType::VectorType;
    using Float = typename VectorType::ValueType;
    using Index = size_t;

    using Mat33 = algebra::Matrix<3,3>;

    /**
     * Default constructor.
     * This will initialize a basis linear hexahedron ( -1 < x < 1, -1 < y < 1, -1 < z < 1) centered
     * on zero (0,0,0) and of dimension 2x2x2.
     */
    LinearRegularHexahedron() : m_anchor_node({-1, -1, -1}), m_H({2, 2, 2})
    {};

    /**
     * Constructor by passing the anchor node (bottom-left-front most node).
     * This will initialize a linear hexahedron of dimension H (hx, hy, hz)
     */
    LinearRegularHexahedron(const PointType & anchor_node, const VectorType & H) : m_anchor_node(anchor_node), m_H(H)
    {};

    /**
     * Constructor by a list initializer of PointType.
     *
     * Example:
     *
     * \code{.cpp}
     * Point3D p1, p2, p3, p4, p5, p6, p7, p8;
     * LinearRegularHexahedron hexahedron ({p1, p2, p3, p4, p5, p6, p7, p8});
     * \endcode
     *
     */
    template <typename OtherPointType>
    LinearRegularHexahedron (const std::initializer_list<OtherPointType> & il)
    {
        // Create a temporary limear (non regular) hexa from the initializer list to compute the hx, hy and hz values
        LinearHexahedron hexa (il);
        m_anchor_node = hexa.node(0);
        m_H = {
                hexa.edge(0).length(), // hx
                hexa.edge(1).length(), // hy
                hexa.edge(8).length()  // hz
        };
    }

    /** Constructor by an array of points **/
    template <typename OtherPointType>
    LinearRegularHexahedron (const std::array<OtherPointType, 8> & n)
    {
        // Create a temporary limear (non regular) hexa from the initializer list to compute the hx, hy and hz values
        LinearHexahedron hexa (n);
        m_anchor_node = hexa.node(0);
        m_H = {
                hexa.edge(0).length(), // hx
                hexa.edge(1).length(), // hy
                hexa.edge(8).length()  // hz
        };
    }

    /** Copy constructor **/
    LinearRegularHexahedron (const LinearHexahedron & other) {
        m_anchor_node = other.node(0);
        m_H = {
                other.edge(0).length(), // hx
                other.edge(1).length(), // hy
                other.edge(8).length()  // hz
        };
    }

    /** Copy constructor **/
    LinearRegularHexahedron (const LinearRegularHexahedron & other)
    : m_anchor_node(other.m_anchor_node)
    , m_H(other.m_H)
    {
    }

    /**
     * Get the node at index
     *
     *
     *          7----------------6
     *         /|               /|
     *        / |              / |
     *       /  |             /  |
     *      /   |            /   |
     *     3----------------2    |
     *     |    |           |    |
     *     |    |           |    |
     *     |    |           |    |
     *     |    4-----------|----5
     *     |   /            |   /
     *     |  /             |  /
     *     | /              | /
     *     |/               |/
     *     0----------------1
     *
     */
    PointType node(Index index) const {
        const auto & hx = m_H[0];
        const auto & hy = m_H[1];
        const auto & hz = m_H[2];

        switch(index) {
            case 0 : return m_anchor_node;
            case 1 : return m_anchor_node + VectorType({hx, 0,  0});
            case 2 : return m_anchor_node + VectorType({hx, hy, 0});
            case 3 : return m_anchor_node + VectorType({0,  hy, 0});
            case 4 : return m_anchor_node + VectorType({0,  0,  hz});
            case 5 : return m_anchor_node + VectorType({hx, 0,  hz});
            case 6 : return m_anchor_node + VectorType({hx, hy, hz});
            case 7 : return m_anchor_node + VectorType({0,  hy, hz});
        }

        throw std::out_of_range("Trying to access the node #" + std::to_string(index) + " of a linear hexahedrons that only contains 8 nodes.");
    };

    /** Return a scaled copy of this hexahedron by s (sx, sy, sz) from the origin **/
    inline LinearRegularHexahedron
    scaled(VectorType s) const
    {
        const PointType anchor = m_anchor_node.scaled(s);
        const VectorType H = m_H.direct_mult(s);

        return LinearRegularHexahedron (anchor, H);
    };

    /** Return a scaled copy of this hexahedron by s from the origin **/
    inline LinearRegularHexahedron
    scaled(Float s) const
    {
        const PointType anchor = m_anchor_node.scaled(s);
        const VectorType H = m_H * s;

        return LinearRegularHexahedron (anchor, H);
    }

    /** Scale this hexahedron by s (sx, sy, sz) from the origin **/
    inline LinearRegularHexahedron &
    scale(VectorType s)
    {
        m_anchor_node.scale(s);
        m_H = m_H.direct_mult(s);

        return (*this);
    }

    /** Scale this hexahedron by s from the origin **/
    inline LinearRegularHexahedron &
    scale(Float s)
    {
        m_anchor_node.scale(s);
        m_H = m_H * s;

        return (*this);
    }

protected:
    PointType m_anchor_node; ///< Anchor node of the hexahedron. This is the first (bottom-left-front) node.
    VectorType m_H; ///< Dimensions (hx, hy and hz) of the hexahedron
};

} // namespace geometry

} // namespace caribou
#endif //CARIBOU_GEOMETRY_LINEARREGULARHEXAHEDRON_H
