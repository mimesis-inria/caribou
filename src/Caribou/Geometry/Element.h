#pragma once

#include <Caribou/config.h>
#include <Caribou/traits.h>
#include <Caribou/macros.h>
#include <Eigen/Dense>
#include <Eigen/QR>
#include <vector>

#include <iostream>

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

template<typename Derived, typename ScalarType = FLOATING_POINT_TYPE>
struct Element {
    // Types
    using Scalar = ScalarType;

    template <INTEGER_TYPE Dim>
    using Vector = Eigen::Matrix<Scalar, Dim, 1>;

    template <INTEGER_TYPE Rows, INTEGER_TYPE Cols, int Options = Eigen::ColMajor>
    using Matrix = Eigen::Matrix<Scalar, Rows, Cols, Options>;
    using Mat3x3   = Eigen::Matrix<double, 3, 3>;

    template <INTEGER_TYPE Rows, INTEGER_TYPE Cols, int Options = Eigen::ColMajor>
    using MatrixI = Eigen::Matrix<Scalar, Rows, Cols, Options>;

    static constexpr auto CanonicalDimension = traits<Derived>::CanonicalDimension;
    static constexpr auto Dimension = traits<Derived>::Dimension;
    static constexpr auto NumberOfNodesAtCompileTime = traits<Derived>::NumberOfNodesAtCompileTime;
    static constexpr auto NumberOfGaussNodesAtCompileTime = traits<Derived>::NumberOfGaussNodesAtCompileTime;
    using LocalCoordinates = Vector<CanonicalDimension>;
    using WorldCoordinates = Vector<Dimension>;

    struct GaussNode {
        LocalCoordinates position;
        Scalar weight;
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
            return 0;
        }
    }

    inline auto local_base() const -> Mat3x3 {
        const auto ex = world_coordinates({1, 0, 0}); 
        const auto ey = world_coordinates({0, 1, 0}); 
        const auto ez = world_coordinates({0, 0, 1});

        Mat3x3 base; 
        base << ex, ey, ez;
        return base; 
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
        static_assert(element_has_boundaries_v<Derived>, "This element type has no boundary elements defined.");
        return self().get_boundary_elements_nodes();
    }

    /**
     * Construct and return the given boundary element.
     *
     * Example:
     * \code{.cpp}
     * Hexahedron<Linear> hexa;
     * Quad<Linear, _3D> face_0 = hexa.boundary_element(0);
     *
     * Tetrahedron<Quadratic> tetra;
     * Triangle<Quadratic, _3D> face_2 = tetra.boundary_element(2);
     * \endcode
     *
     */
    inline auto boundary_element(const UNSIGNED_INTEGER_TYPE & boundary_id) const {
        static_assert(element_has_boundaries_v<Derived>, "This element type has no boundary elements defined.");
        using BoundaryElement = typename traits<Derived>::BoundaryElementType;
        const auto & node_indices = boundary_elements_node_indices()[boundary_id];
        Matrix<BoundaryElement::NumberOfNodesAtCompileTime, Dimension> nodes;
        for (Eigen::Index boundary_node_id = 0; boundary_node_id < nodes.rows(); ++boundary_node_id) {
            nodes.row(boundary_node_id) = self().node(node_indices[boundary_node_id]);
        }
        return BoundaryElement(nodes);
    }

    /**
     * Get the Lagrange polynomial values evaluated at local coordinates xi w.r.t each element's interpolation nodes.
     *
     * Example:
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
     * Example:
     * \code{.cpp}
     * // Computes the derivatives of node #2 Lagrange polynomial evaluated at local coordinates {-0.4}
     * float dp = Segment2::dL(-0.4)[2];
     * \endcode
     */
    inline auto dL(const LocalCoordinates & xi) const -> Matrix<NumberOfNodesAtCompileTime, CanonicalDimension> {return self().get_dL(xi);};

    /** Get the position at the center of the element */
    inline auto center() const -> WorldCoordinates {return self().get_center();}

    /** Get the world coordinates of a point from its local coordinates. */
    inline auto world_coordinates(const LocalCoordinates & coordinates) const -> WorldCoordinates {
        return WorldCoordinates(self().interpolate(coordinates, self().nodes()));
    }

    /**
     * Get the local coordinates of a point from its world coordinates by doing a set of Newton-Raphson iterations.
     *
     * \sa local_coordinates() for more details.
     *
     * \note By default, the Newton-Raphson will start by an approximation of the local coordinates at [0, 0, 0].
     *       The iterations will stop at 5 iterations, or if the norm o*       f relative residual |R|/|R0| is less than 1e-5.
     */
    inline auto local_coordinates(const WorldCoordinates & coordinates) const -> LocalCoordinates {
        return local_coordinates(coordinates, LocalCoordinates::Constant(0), 1e-5, 5);
    }

    /**
     * Get the local coordinates of a point from its world coordinates by doing a set of Newton-Raphson iterations.
     *
     * By taking the Taylor expansion of the transformation \f$T(\vec{\xi}) \rightarrow \vec{x} \f$
     * with \f$\vec{x} = [x,y,z]^T\f$ being the world coordinates of a point and \f$\vec{\xi} = [u,v,w]^T\f$ its local coordinates
     * within the element, we have
     *
     * \f{eqnarray*}{
     *     x_p &= x_0 + \frac{\partial x}{\partial u} \cdot (u_p - u_0) + \frac{\partial x}{\partial v} \cdot (v_p - v_0) + \frac{\partial x}{\partial w} \cdot (w_p - w_0) \\
     *     y_p &= y_0 + \frac{\partial y}{\partial u} \cdot (u_p - u_0) + \frac{\partial y}{\partial v} \cdot (v_p - v_0) + \frac{\partial y}{\partial w} \cdot (w_p - w_0) \\
     *     z_p &= z_0 + \frac{\partial z}{\partial u} \cdot (u_p - u_0) + \frac{\partial z}{\partial v} \cdot (v_p - v_0) + \frac{\partial z}{\partial w} \cdot (w_p - w_0) \\
     * \f}
     *
     * where partial derivatives are evaluated at \f$\vec{\xi}_0 = [u_0,v_0,w_0]^T\f$. We can reformulate with the following matrix form
     *
     * \f{eqnarray*}{
     *     \vec{x}_p = \vec{x}_0 + \mathrm{J} \cdot (\vec{\xi}_p - \vec{\xi}_0)
     * \f}
     *
     * where \f$ \vec{x}_p = T(\vec{\xi}_p) \f$, \f$ \vec{x}_0 = T(\vec{\xi}_0) \f$ and \f$ \mathrm{J} \f$ is the Jacobian
     * of the transformation \f$T\f$. Since we are trying to find \f$\vec{\xi}_p\f$, we can rearange the last equation to get
     *
     * \f{eqnarray*}{
     *     \vec{\xi}_p = \vec{\xi}_0 + \mathrm{J}^{-1} (\vec{x}_p - \vec{x}_0)
     * \f}
     *
     * Hence, starting from an initial guess at local coordinates \f$ \vec{\xi}_0 \f$, we have the following iterative
     * method:
     *
     * \f{eqnarray*}{
     *     \vec{\xi}_p^{k+1} = \vec{\xi}_k + \mathrm{J}^{-1} (\vec{x}_p - T(\vec{\xi}_k))
     * \f}
     *
     * The iterations stop when \f$ \frac{|\vec{x}_p - T(\vec{\xi}_k)|}{|\vec{x}_p - T(\vec{\xi}_0)|} < \epsilon \f$
     *
     * \note
     * When trying to find the local coordinates of non-matching manifolds (for example, the local coordinates
     * of a triangle in a 3D manifold), the following recursive formulae is used:
     * \f{eqnarray*}{
     *     \vec{\xi}_p^{k+1} = \vec{\xi}_k + (\mathrm{J}^T\mathrm{J})^{-1} \mathrm{J}^T (\vec{x}_p - T(\vec{\xi}_k))
     * \f}
     *
     * @param coordinates    The world coordinates of the point from which we want to get the local coordinates.
     * @param starting_point An approximation of the real local coordinates we want the get. The closer it is to the
     *                       solution, the faster the Newton will converge.
     * @param residual_tolerance The threshold of relative norm of the residual at which point the Newton is
     *                           said to converge (|R|/|R0| < threshold).
     * @param maximum_number_of_iterations The maximum number of Newton-Raphson iterations we can do before divergence.
     * @return The local coordinates of the point at the last Newton-Raphson iteration completed.
     */
    inline auto local_coordinates(
        const WorldCoordinates & coordinates,
        const LocalCoordinates & starting_point,
        const FLOATING_POINT_TYPE & residual_tolerance,
        const UNSIGNED_INTEGER_TYPE & maximum_number_of_iterations) const -> LocalCoordinates {

        using namespace Eigen;

        LocalCoordinates xi = starting_point;
        UNSIGNED_INTEGER_TYPE iteration = 0;
        FLOATING_POINT_TYPE squared_threshold = residual_tolerance*residual_tolerance;

        // Initial residual
        WorldCoordinates residual = coordinates - world_coordinates(xi);
        FLOATING_POINT_TYPE r_norm_2 = residual.squaredNorm();
        FLOATING_POINT_TYPE r0_norm_2 = r_norm_2;

        if (r_norm_2 < squared_threshold) {
            // The initial guess is good enough
            return xi;
        }

        // Start the iterations
        do {
            Matrix<Dimension, CanonicalDimension> J = jacobian(xi);

            LocalCoordinates dxi;
            if constexpr (Dimension == CanonicalDimension) {
                dxi.noalias() = J.inverse() * residual;
            } else {
                dxi.noalias() = (J.transpose()*J).inverse() * (J.transpose()*residual);
            }

            xi.noalias() = (xi + dxi).eval();
            residual.noalias() = coordinates - world_coordinates(xi);
            r_norm_2 = residual.squaredNorm();

            ++iteration;
        }
        while (iteration < maximum_number_of_iterations and r_norm_2 >= r0_norm_2*squared_threshold);

        return xi;
    }

    /**
     * Return true if the element contains the point located at the given local coordinates.
     * @param coordinates Local coordinates of a point
     * @param eps If the given point is located barely outside the element, which is, less than this eps value, returns true.
     */
    inline auto contains_local(const LocalCoordinates & xi, const FLOATING_POINT_TYPE & eps = 1e-10) const -> bool {
        return self().get_contains_local(xi, eps);
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
        using MatrixScalar = typename Eigen::MatrixBase<MatrixType>::Scalar;
        const auto result = ((values.array().colwise() * self().L(coordinates).array().template cast<typename Eigen::MatrixBase<MatrixType>::Scalar>()).matrix().colwise().sum().transpose()).eval();
        if constexpr (NbCols == 1) {
            return static_cast<MatrixScalar>(result[0]);
        } else {
            return result.template cast<MatrixScalar>();
        }
    }

    /**
     * Compute the Jacobian matrix of the transformation T(xi)-> x evaluated at local coordinates xi.
     *
     * The Jacobian is defined as:
     *
     * \verbatim
     *
     * |dx|     |du|
     * |dy| = J |dv|
     * |dz|     |dw|
     *
     * 1D canonical element
     * --------------------
     *
     * 1D manifold:    J(u) = dx/du = sum_i dNi/du * x_i
     *                 det(J) = |J|
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
     *                 det(J) = |det(J)|
     *
     *                          | dx/du   dx/dv |   | sum dNi/du  x_i    sum dNi/dv  x_i |
     * 3D manifold:    J(u,v) = | dy/du   dy/dv | = | sum dNi/du  y_i    sum dNi/dv  y_i |
     *                          | dz/du   dz/dv | = | sum dNi/du  z_i    sum dNi/dv  z_i |
     *                 det(J) = sqrt((J.transpose() * J).determinant())
     *                        = J.col(0).cross(J.col(1)).norm();
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
     *\endverbatim
     *
     * Example:
     * \code{.cpp}
     * // Computes the Jacobian of a 3D segment and its determinant evaluated at local coordinates 0.5 (half-way through the segment)
     * Segment<3, Linear> segment {{5, 5, 5}, {10, 5,0}};
     * Matrix<3,1> J = segment.jacobian (0.5);
     * double detJ = sqrt(J.dot(J));
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