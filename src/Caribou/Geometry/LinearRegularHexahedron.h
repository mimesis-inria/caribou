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

struct LinearRegularHexahedron
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
        static constexpr std::array<std::array<unsigned char, 2>, 12> edges
        {{
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
                node(edges[index][0]),
                node(edges[index][1])
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
        static const std::array<std::array<unsigned char, 4>, 6> nodes_of_face =
        {{
            {0, 1, 2, 3}, // 0 ; Edges = {0, 1, 2, 3,}
            {3, 7, 4, 0}, // 1 ; Edges = {11, 7, 8, 3}
            {7, 6, 5, 4}, // 2 ; Edges = {6, 5, 4, 7}
            {5, 6, 2, 1}, // 3 ; Edges = {5, 10, 1, 9}
            {5, 1, 0, 4}, // 4 ; Edges = {9, 0, 8, 4}
            {2, 6, 7, 3}, // 5 ; Edges = {10, 6, 11, 2}
        }};

        return FaceType {
                node(nodes_of_face[index][0]),
                node(nodes_of_face[index][1]),
                node(nodes_of_face[index][2]),
                node(nodes_of_face[index][3])
        };
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
    };

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
    };

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
    };

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
        Mat33 result (true /* initialize_to_zero */);

        for (Index i = 0; i < 8; ++i) {
            VectorType coordinates = node(i);
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
    };

    /**
     * Compute an integral approximation by gauss quadrature on the hexahedron of the given evaluation function.
     *
     * @example
     * \code{.cpp}
     * // Integrate the polynomial 1 + 2x + 2xy + 3*z on an hexahedron.
     * float result = LinearHexahedron(x1, x2, x3, x4, x5, x6, x7, x8).gauss_integrate(
     *   [] (const LinearHexahedron & hexa, const float & xi, const float & eta, const float & zeta) -> float {
     *     return 1 + 2*xi + 2*xi*eta + 3*zeta;
     *   }
     * );
     * \endcode
     *
     * @tparam ValueType The result type of the evaluation function.
     * This type must implement the assignment "=", assignment-addition "+=", and multiplication "*" with a scalar type (float, double) operators.
     * @tparam EvaluateFunctionType Callback function reference type. See evaluate parameter.
     *
     * @param evaluate
     * Callback function of the signature
     *
     *     ValueType f (const LinearHexahedron & hexa, const float & xi, const float & eta, const float & zeta);
     *
     * Where hexa is a reference to the current hexahadron on which we integrate, and the coordinates xi, eta and zeta
     * forms the local position of a sample point on which we want to get the evaluation value of type ValueType.
     *
     * @return The value of the integral computed on this hexahedron.
     *
     */
    template <typename ValueType , typename EvaluateFunctor>
    ValueType gauss_quadrature(const ValueType & initial_value, EvaluateFunctor evaluate) const
    {
        const static std::array<VectorType, 8> gauss_points =
        {{
            {-1./sqrt(3.0), -1./sqrt(3.0), -1./sqrt(3.0)},
            {-1./sqrt(3.0), -1./sqrt(3.0),  1./sqrt(3.0)},
            {-1./sqrt(3.0),  1./sqrt(3.0), -1./sqrt(3.0)},
            {-1./sqrt(3.0),  1./sqrt(3.0),  1./sqrt(3.0)},
            { 1./sqrt(3.0), -1./sqrt(3.0), -1./sqrt(3.0)},
            { 1./sqrt(3.0), -1./sqrt(3.0),  1./sqrt(3.0)},
            { 1./sqrt(3.0),  1./sqrt(3.0), -1./sqrt(3.0)},
            { 1./sqrt(3.0),  1./sqrt(3.0),  1./sqrt(3.0)}
         }};

        ValueType result = initial_value;

        for (const VectorType p : gauss_points) {
            const Float & xi   = p[0];
            const Float & eta  = p[1];
            const Float & zeta = p[2];

            Float detJ = Jacobian(xi, eta, zeta).determinant();
            result += evaluate(*this, xi, eta, zeta) * detJ;
        }

        return result;
    }

    /** Get the global (world) coordinates (x, y, z) of a point inside the hexa from its local coordinates (xi, eta, zeta). */
    inline VectorType
    from_local_coordinate(const VectorType & local_coordinate) {
        const Float & xi =   local_coordinate[0];
        const Float & eta =  local_coordinate[1];
        const Float & zeta = local_coordinate[2];

        return (
                N(0, xi, eta, zeta) * node(0) +
                N(1, xi, eta, zeta) * node(1) +
                N(2, xi, eta, zeta) * node(2) +
                N(3, xi, eta, zeta) * node(3) +
                N(4, xi, eta, zeta) * node(4) +
                N(5, xi, eta, zeta) * node(5) +
                N(6, xi, eta, zeta) * node(6) +
                N(7, xi, eta, zeta) * node(7)
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

protected:
    PointType m_anchor_node; ///< Anchor node of the hexahedron. This is the first (bottom-left-front) node.
    VectorType m_H; ///< Dimensions (hx, hy and hz) of the hexahedron
};

} // namespace geometry

} // namespace caribou
#endif //CARIBOU_GEOMETRY_LINEARREGULARHEXAHEDRON_H
