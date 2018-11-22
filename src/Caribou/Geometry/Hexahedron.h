#ifndef CARIBOU_GEOMETRY_HEXAHEDRON_H
#define CARIBOU_GEOMETRY_HEXAHEDRON_H

#include <Caribou/Geometry/Point.h>
#include <Caribou/Geometry/Segment.h>
#include <Caribou/Geometry/Quad.h>

#include <Caribou/Algebra/Matrix.h>

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
template <size_t NNodes>
struct Hexahedron
{

    static constexpr size_t NumberOfNodes = NNodes;

    using PointType = Point<3>;
    using VectorType = typename PointType::VectorType;
    using ValueType = typename VectorType::ValueType;
    using SegmentType = Segment<3>;
    using Mat33 = algebra::Matrix<3,3>;

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

    /** Constructor by an array of points **/
    Hexahedron (const std::array<PointType, NumberOfNodes> & n) : nodes(n)
    {
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

    /** Copy constructor **/
    explicit Hexahedron (const Hexahedron<NumberOfNodes> & other) {
        std::copy(std::begin(other.nodes), std::end(other.nodes), std::begin(nodes));
    }

    /** Scale the hexahedron by s (sx, sy, sz) from the origin **/
    Hexahedron<NumberOfNodes>
    inline scale(VectorType s) const
    {
        std::array<PointType, NumberOfNodes> scaled_nodes;
        for (size_t i = 0; i < NumberOfNodes; ++i) {
            scaled_nodes[i] = nodes[i].scale(s);
        }

        return scaled_nodes;
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

struct LinearHexahedron : public Hexahedron<8>
{
    using Base = Hexahedron<8>;
    using Hexahedron<8>::Hexahedron; // Import the base constructors
    using VectorType = typename Base::VectorType;
    using PointType = typename Base::PointType;
    using SegmentType = typename Base::SegmentType;
    using FaceType = Quad<3>;
    using Float = typename VectorType::ValueType;
    using Index = size_t;

    LinearHexahedron() : Base ({
       //  {xi, eta, zeta}
       PointType({-1, -1, -1}), // 0
       PointType({ 1, -1, -1}), // 1
       PointType({ 1,  1, -1}), // 2
       PointType({-1,  1, -1}), // 3
       PointType({-1, -1,  1}), // 4
       PointType({ 1, -1,  1}), // 5
       PointType({ 1,  1,  1}), // 6
       PointType({-1,  1,  1})  // 7
    })
    {}

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
        return Base::nodes[index];
    };

    /**
     * Get the edge at index
     * 
     *      7            6
     *       .           .
     *        . +----------------+
     *         /|               /|
     * 11 ..../ |              / |
     *       /  |  2          /....... 10
     *      /   |  .         /   |
     *     +----------------+    |.... 5
     *     |    |           |    |
     *     |    |           |    |
     *     |    |           |......... 1
     * 3 ..|    +-----------|----+
     *     |   /      4     |   /
     * 8 ..|../             |  /...... 9
     *     | /              | /
     *     |/        0      |/
     *     +----------------+
     *
     */
    SegmentType edge(Index index) const {
        static constexpr std::array<std::array<unsigned char, 2>, 12> edges {{
            {0, 1}, // 0
            {1, 2}, // 1
            {2, 3}, // 2
            {3, 0}, // 3
            {4, 5}, // 4
            {5, 6}, // 5
            {6, 7}, // 6
            {7, 4}, // 7
            {0, 4}, // 8
            {1, 5}, // 9
            {2, 6}, // 10
            {3, 7}, // 11
        }};

        return SegmentType (
                Base::nodes[edges[index][0]],
                Base::nodes[edges[index][1]]
        );
    };

    /**
     * Get the face at index
     *                          2
     *                         .
     *          +----------------+
     *         /|     5         /|
     *        / |     .        / |
     *       /  |     .       /  |
     *      /   |            /   |
     *     +----------------+    |
     *     |    |           |    |
     *     |    |           |  .......... 3
     * 1...|    |           |    |
     *     |    +-----------|----+
     *     |   /      ................... 0
     *     |  /             |  /
     *     | /              | /
     *     |/               |/
     *     +----------------+
     *               .
     *               .
     *               4
     */
    FaceType face(Index index) const {
        static const std::array<std::array<unsigned char, 4>, 6> nodes_of_face = {{
            {0, 1, 2, 3}, // 0 ; Edges = {0, 1, 2, 3,}
            {3, 7, 4, 0}, // 1 ; Edges = {11, 7, 8, 3}
            {7, 6, 5, 4}, // 2 ; Edges = {6, 5, 4, 7}
            {5, 6, 2, 1}, // 3 ; Edges = {5, 10, 1, 9}
            {5, 1, 0, 4}, // 4 ; Edges = {9, 0, 8, 4}
            {2, 6, 7, 3}, // 5 ; Edges = {10, 6, 11, 2}
        }};

        return FaceType {
                Base::nodes[nodes_of_face[index][0]],
                Base::nodes[nodes_of_face[index][1]],
                Base::nodes[nodes_of_face[index][2]],
                Base::nodes[nodes_of_face[index][3]]
         };
    };

    /**
     * Compute the shape value of the hexa's node i evaluated at local coordinates {xi, eta, zeta}.
     */
    inline Float
    N(const Index i, const Float & xi, const Float & eta, const Float & zeta) const
    {
        const VectorType Xi_i = LinearHexahedron().node(i); // Local coordinates of each node
        const Float & xi_i =   Xi_i[0];
        const Float & eta_i =  Xi_i[1];
        const Float & zeta_i = Xi_i[2];

        return Float(1/8.) * (1 + xi_i*xi) * (1 + eta_i*eta) * (1 + zeta_i*zeta);
    }

    /**
     * Compute the shape (N) derivatives w.r.t the local frame Xi {dN/dxi, dN/deta, dN/dzeta} of the hexa's node i evaluated
     * at local coordinates {xi, eta, zeta}.
     */
    inline VectorType
    dN_dXi(const Index i, const Float & xi, const Float & eta, const Float & zeta) const
    {
        const VectorType Xi_i = LinearHexahedron().node(i); // Local coordinates of each node
        const Float & xi_i =   Xi_i[0];
        const Float & eta_i =  Xi_i[1];
        const Float & zeta_i = Xi_i[2];

        return VectorType({
           1./8 * (        xi_i     * (1 + eta_i*eta) * (1 + zeta_i*zeta)   ), // dN / dxi
           1./8 * (   (1 + xi_i*xi) *      eta_i      * (1 + zeta_i*zeta)   ), // dN / deta
           1./8 * (   (1 + xi_i*xi) * (1 + eta_i*eta) *      zeta_i         )  // dN / dzeta
        });
    }

    /**
     * Compute the shape (N) derivatives w.r.t. the global frame X (dN/dx, dN/dy, dN/dz) of the hexa's node i evaluated
     * at local coordinates {xi, eta, zeta}.
     */
    inline VectorType
    dN_dX(const Index i, const Float & xi, const Float & eta, const Float & zeta) const
    {
        Mat33 J = Jacobian(xi, eta, zeta);
        const VectorType dNdXi = dN_dXi(i, xi, eta, zeta);

        return (J^-1) * dNdXi;
    }

    /**
     * Compute the Jacobian matrix of the shape function N evaluated at local coordinates {xi, eta, zeta}.
     *
     * The Jacobian is defined as:
     *
     *     | dx/dXi    dy/dXi   dz/dXi   |   | sum dNi/dXi   xi    sum dNi/dXi   yi    sum dNi/dXi   zi |
     * J = | dx/dEta   dy/dEta  dz/dEta  | = | sum dNi/dEta  xi    sum dNi/dEta  yi    sum dNi/dEta  zi |
     *     | dx/dZeta  dy/dZeta dz/dZeta |   | sum dNi/dZeta xi    sum dNi/dZeta yi    sum dNi/dZeta zi |
     *
     * where dN/dxi (resp. deta, dzeta) is the partial derivative of the shape function at node i (xi, yi, zi)
     * w.r.t the local frame Xi {Xi, Eta, Zeta} evaluated at local coordinate  {xi, eta, zeta}
     */
    inline Mat33
    Jacobian(const Float & xi, const Float & eta, const Float & zeta) const
    {
        Mat33 result (true /* initizalise_to_zero */);

        for (Index i = 0; i < NumberOfNodes; ++i) {
            VectorType coordinates = node(i).coordinates;
            const Float & x = coordinates[0];
            const Float & y = coordinates[1];
            const Float & z = coordinates[2];
            const VectorType dNi_dXi = dN_dXi(i, xi, eta, zeta);
            Mat33 j = {{
                {dNi_dXi[0]*x, dNi_dXi[0]*y, dNi_dXi[0]*z},
                {dNi_dXi[1]*x, dNi_dXi[1]*y, dNi_dXi[1]*z},
                {dNi_dXi[2]*x, dNi_dXi[2]*y, dNi_dXi[2]*z}
            }};
            result += j;
        };

        return result;
    }

    /** Get the global (world) coordinates (x, y, z) of a point inside the hexa from its local coordinates (xi, eta, zeta). */
    inline VectorType
    from_local_coordinate(const VectorType & local_coordinate) {
        const Float & xi =   local_coordinate[0];
        const Float & eta =  local_coordinate[1];
        const Float & zeta = local_coordinate[2];
        std::array<VectorType, 8> x = {
                node(0).coordinates, node(1).coordinates, node(2).coordinates, node(3).coordinates,
                node(4).coordinates, node(5).coordinates, node(6).coordinates, node(7).coordinates
        };

        return (
            N(0, xi, eta, zeta) * x[0] +
            N(1, xi, eta, zeta) * x[1] +
            N(2, xi, eta, zeta) * x[2] +
            N(3, xi, eta, zeta) * x[3] +
            N(4, xi, eta, zeta) * x[4] +
            N(5, xi, eta, zeta) * x[5] +
            N(6, xi, eta, zeta) * x[6] +
            N(7, xi, eta, zeta) * x[7]
        );
    };

    inline FaceType
    left() const {return face(1);};

    inline FaceType
    right() const {return face(3);};

    inline FaceType
    top() const {return face(5);};

    inline FaceType
    bottom() const {return face(4);};

    inline FaceType
    back() const {return face(2);};

    inline FaceType
    front() const {return face(0);};
};

/**
 * Regular linear hexahedron
 * @tparam TVector
 */
struct RegularLinearHexahedron : public LinearHexahedron
{
    using LinearHexahedron::LinearHexahedron; // Import the base constructors
};

} // namespace geometry

} // namespace caribou

#endif //CARIBOU_GEOMETRY_HEXAHEDRON_H
