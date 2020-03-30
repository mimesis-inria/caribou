#pragma once

#include <Caribou/config.h>
#include <Caribou/Traits.h>
#include <Caribou/macros.h>
#include <Eigen/Core>

namespace caribou::geometry {

template<typename T>
struct traits;

/**
 * Test whether or not the Element type T has boundary elements.
 *
 * This is done by looking into the traits<T> and verify if it as a type named BoundaryElementType.
 */
template<class T, class Enable = void>
struct element_has_boundaries : std::false_type {};

template< class T>
constexpr bool element_has_boundaries_v = element_has_boundaries<T>::value;

template<typename Derived>
struct Element {
    // Types
    template <UNSIGNED_INTEGER_TYPE Dim>
    using Vector = Eigen::Matrix<FLOATING_POINT_TYPE, Dim, 1>;

    template <UNSIGNED_INTEGER_TYPE Rows, UNSIGNED_INTEGER_TYPE Cols, int Options = 0>
    using Matrix = Eigen::Matrix<FLOATING_POINT_TYPE, Rows, Cols, Options>;

    static constexpr auto CanonicalDimension = traits<Derived>::CanonicalDimension;
    static constexpr auto Dimension = traits<Derived>::Dimension;
    static constexpr auto NumberOfNodesAtCompileTime = traits<Derived>::NumberOfNodesAtCompileTime;
    static constexpr auto NumberOfGaussNodesAtCompileTime = traits<Derived>::NumberOfGaussNodesAtCompileTime;
    using LocalCoordinates = Vector<CanonicalDimension>;
    using WorldCoordinates = Vector<Dimension>;

    struct GaussNode {
        LocalCoordinates position;
        FLOATING_POINT_TYPE weight;
    };


    // Functions
    /** Get the number of nodes in the element */
    [[nodiscard]]
    inline auto number_of_nodes() const -> UNSIGNED_INTEGER_TYPE {return self().get_number_of_nodes();}

    /** Get the number of gauss nodes in the element */
    [[nodiscard]]
    inline auto number_of_gauss_nodes() const -> UNSIGNED_INTEGER_TYPE {return self().get_number_of_gauss_nodes();}

    /** Get the Node at given index */
    inline auto node(const UNSIGNED_INTEGER_TYPE & index) const {return self().get_node(index);};

    /** Get the set of nodes */
    inline auto nodes() const {return self().get_nodes();};

    /** Get the gauss node at given index */
    inline auto gauss_node(const UNSIGNED_INTEGER_TYPE & index) const -> const GaussNode & {return self().gauss_nodes()[index];};

    /** Get the set of gauss nodes */
    inline auto gauss_nodes() const -> const std::vector<GaussNode> & {return self().get_gauss_nodes();};

    /** Get the number of boundary elements (ex. number of triangles in a tetrahedron) */
    inline auto number_of_boundary_elements() const {
        if constexpr (element_has_boundaries_v<Derived>) {
            return self().get_number_of_boundary_elements();
        } else {
            static_assert(element_has_boundaries_v<Derived>, "This element type has no boundary elements defined.");
        }
    }

    /**
     * Get the list of node indices of the boundary elements.
     *
     * The return type is up to the implementation but should normally be:
     *   1. Dynamic number of boundary element where each boundary element has a dynamic number of nodes
     *      std::vector<std::vector<unsigned int>>
     *   2. Dynamic number of boundary element where each boundary element has a static number of nodes
     *      std::vector<std::array<unsigned int, traits<BoundaryElement>::NumberOfNodesAtCompileTime>>
     *   3. Static number of boundary element where each boundary element has a dynamic number of nodes
     *      std::array<std::vector<unsigned int>, traits<Element>::NumberOfBoundaryElementsAtCompileTime>
     *   4. Static number of boundary element where each boundary element has a static number of nodes
     *   std::array<
     *     std::array<unsigned int, traits<BoundaryElement>::NumberOfNodesAtCompileTime>,
     *     traits<Element>::NumberOfBoundaryElementsAtCompileTime
     *   >
     */
    inline auto boundary_elements_node_indices() const -> const auto & {
        if constexpr (element_has_boundaries_v<Derived>) {
            return self().get_boundary_elements_nodes();
        } else {
            static_assert(element_has_boundaries_v<Derived>, "This element type has no boundary elements defined.");
        }
    }

    /**
     * Get the Lagrange polynomial values evaluated at local coordinates xi w.r.t each element's interpolation nodes.
     *
     * @example
     * \code{.cpp}
     * // Computes the value of node #2 Lagrange polynomial evaluated at local coordinates {-0.4} on a segment
     * Segment<3, Linear> segment;
     * float p = segment.L(-0.4)[2];
     * \endcode
     */
    inline auto  L(const LocalCoordinates & xi) const -> Vector<NumberOfNodesAtCompileTime> {return self().get_L(xi);};

    /**
     * Get the Lagrange polynomial derivatives w.r.t the local frame {dL/du} evaluated at local
     * coordinates {u} w.r.t each segment's interpolation nodes.
     *
     * @example
     * \code{.cpp}
     * // Computes the derivatives of node #2 Lagrange polynomial evaluated at local coordinates {-0.4}
     * float dp = Segment2::dL(-0.4)[2];
     * \endcode
     */
    inline auto dL(const LocalCoordinates & xi) const -> Matrix<NumberOfNodesAtCompileTime, CanonicalDimension> {return self().get_dL(xi);};

    /** Get the position at the center of the element */
    inline auto center() const -> WorldCoordinates {return self().get_center();}

    /** Get the world coordinates of a point from its local coordinates. */
    inline auto world_coordinates(const LocalCoordinates & coordinates) const {
        return WorldCoordinates(self().interpolate(coordinates, self().nodes()));
    }

    /**
     * Interpolate a value at local coordinates from the given interpolation node values.
     *
     * The values at nodes must be provided in an Eigen::Matrix where the row i of the matrix is the value (as a row-vector or a scalar)
     * at the node i.
     *
     * @tparam MatrixType Type of the Eigen matrix containing the values at nodes.
     */
    template <typename MatrixType>
    inline auto interpolate (const LocalCoordinates & coordinates,
                             const Eigen::MatrixBase<MatrixType> & values) const {
        static_assert(Eigen::MatrixBase<MatrixType>::RowsAtCompileTime == NumberOfNodesAtCompileTime,
                      "The matrix containing the values at each nodes must have one node-value per row.");
        constexpr auto NbCols = Eigen::MatrixBase<MatrixType>::ColsAtCompileTime;
        using Scalar = typename Eigen::MatrixBase<MatrixType>::Scalar;
        const auto result = ((values.array().colwise() * self().L(coordinates).array()).matrix().colwise().sum().transpose()).eval();
        if constexpr (NbCols == 1) {
            return static_cast<Scalar>(result[0]);
        } else {
            return result.template cast<Scalar>();
        }
    }

    /**
     * Compute the Jacobian matrix evaluated at local coordinates.
     *
     * The Jacobian is defined as:
     *
     * 1D canonical element
     * --------------------
     *
     * 1D manifold:    J(u) = dx/du = sum_i dNi/du * x_i
     *                 det(J) = J
     *
     * 2D manifold:    J(u)  = | dx/du | = | sum dNi/du x_i |
     *                         | dy/du | = | sum dNi/du y_i |
     *                 det(J) = sqrt(J.dot(J))
     *
     *                        | dx/du | = | sum dNi/du x_i |
     * 3D manifold:    J(u) = | dy/du | = | sum dNi/du y_i |
     *                        | dz/du | = | sum dNi/du z_i |
     *                 det(J) = sqrt(J.dot(J))
     *
     *
     * 2D canonical element
     * --------------------
     *
     * 1D manifold:    N/A
     *
     * 2D manifold:    J(u,v) = | dx/du   dx/dv |   | sum dNi/du  x_i    sum dNi/dv  x_i |
     *                          | dy/du   dy/dv | = | sum dNi/du  y_i    sum dNi/dv  y_i |
     *                 det(J) = det(J)
     *
     *                          | dx/du   dx/dv |   | sum dNi/du  x_i    sum dNi/dv  x_i |
     * 3D manifold:    J(u,v) = | dy/du   dy/dv | = | sum dNi/du  y_i    sum dNi/dv  y_i |
     *                          | dz/du   dz/dv | = | sum dNi/du  z_i    sum dNi/dv  z_i |
     *
     *
     * 3D canonical element
     * --------------------
     *
     * 1D manifold:    N/A
     * 2D manifold:    N/A
     *                            | dx/du   dx/dv   dx/dw |   | sum dNi/du x_i   sum dNi/dv x_i    sum dNi/dw x_i |
     * 3D manifold:    J(u,v,w) = | dy/du   dy/dv   dy/dw | = | sum dNi/du y_i   sum dNi/dv y_i    sum dNi/dw y_i |
     *                            | dz/du   dz/dv   dz/dw |   | sum dNi/du z_i   sum dNi/dv z_i    sum dNi/dw z_i |
     *
     *
     *
     *
     * where dNi/du (resp. dv and dw) is the partial derivative of the shape function at node i
     * w.r.t the local frame of the canonical element evaluated at local coordinate  {u, v, w} and
     * where {xi, yi and zi} are the world coordinates of the position of node i on its element manifold.
     *
     * @example
     * \code{.cpp}
     * // Computes the Jacobian of a 3D segment and its determinant evaluated at local coordinates 0.5 (half-way through the segment)
     * Segment<3, Linear> segment {{5, 5, 5}, {10, 5,0}};
     * Matrix<3,1> J = segment.jacobian (0.5);
     * double detJ = (J^T * J).determinant() ^ 1/2;
     *
     * // Computes the Jacobian of a 3D triangle and its determinant evaluated at local coordinates {1/3, 1/3} (on its center point)
     * Triangle<3, Linear> triangle {{5,5,5}, {15, 5, 5}, {10, 10, 10}};
     * Matrix<3,2> J = triangle.jacobian(1/3., 1/3.);
     * double detJ = (J^T * J).determinant() ^ 1/2;
     *
     * // Computes the Jacobian of a 3D tetrahedron and its determinant evaluated at local coordinates {1/3, 1/3, 1/3} (on its center point)
     * Tetrahedron<3> tetra {{5,5,5}, {15, 5, 5}, {10, 10, 10}};
     * Matrix<3,3> J = tetra.jacobian(1/3., 1/3., 1/3.);
     * double detJ = J.determinant();
     * \endcode
     */
    inline auto jacobian (const LocalCoordinates & coordinates) const -> Matrix<Dimension, CanonicalDimension>
    {
        const auto shape_derivatives = self().dL(coordinates);
        return self().nodes().transpose() * shape_derivatives;
    }
private:
    auto self() -> Derived& { return *static_cast<Derived*>(this); }
    auto self() const -> const Derived& { return *static_cast<const Derived*>(this); }
};


template <class T>
using element_boundary_type_t = typename traits<T>::BoundaryElementType;

template<class T>
struct element_has_boundaries<T, CLASS_REQUIRES(
    caribou::internal::is_detected_v<element_boundary_type_t, T>
)> : std::true_type {};


}  /// namespace caribou::geometry