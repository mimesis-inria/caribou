#pragma once

#include <Caribou/config.h>
#include <Caribou/constants.h>
#include <Caribou/Geometry/Element.h>
#include <Eigen/Core>

namespace caribou::geometry {

template<typename Derived>
struct Beam : public Element<Derived> {
    // Types 
    using Base = Element<Derived>;

    using LocalCoordinates = typename Base::LocalCoordinates;
    using WorldCoordinates = typename Base::WorldCoordinates;

    using GaussNode = typename Base::GaussNode;

    template <UNSIGNED_INTEGER_TYPE Dim>
    using Vector = typename Base::template Vector<Dim>;

    template <UNSIGNED_INTEGER_TYPE Rows, UNSIGNED_INTEGER_TYPE Cols>
    using Matrix = typename Base::template Matrix<Rows, Cols>;

    // Constants
    static constexpr auto CanonicalDimension = 3;
    static constexpr auto Dimension = 3;
    static constexpr auto NumberOfNodesAtCompileTime = 2;
    static constexpr auto NumberOfGaussNodesAtCompileTime = 1;

    using Beam_matrix = typename Matrix<NumberOfNodesAtCompileTime, 3);

    /** Default empty constructor */
    Beam() = default:

    // Constructor with positions and rotations
    Beam(const Beam_matrix positions_matrix, const Beam_matrix rotation_matrix)  {

        // construction of the nodes
        this->p_nodes.row(0) = positions_matrix.row(0);
        this->p_nodes.row(1) = positions_matrix.row(1):

        // construction of the rotations
        this->rotations.row(0) = matrix_rotation.row(0);
        this->rotations.row(1) = matrix_rotation.row(1);
    };

    // Constructor with rotations and defaults positions
    Beam(const Beam_matrix rotation_matrix)  {

        // construction of the nodes
        this->p_nodes.row(0) = WorldCoordinates(+1, 0, 0);
        this->p_nodes.row(1) = WorldCoordinates(-1, 0, 0);

        // construction of the rotations
        this->rotations.row(0) = matrix_rotation.row(0);
        this->rotations.row(1) = matrix_rotation.row(1);
    };

private:
    // Implementations
    friend struct Element<Derived>
    [[nodiscard]]
    inline auto get_number_of_nodes() const {return NumberOfNodesAtCompileTime;}
    inline auto get_number_of_gauss_nodes() const {return NumberOfGaussNodesAtCompileTime;}
    inline auto get_node(const UNSIGNED_INTEGER_TYPE & index) const {return WorldCoordinates(p_nodes.row(index));};
    inline auto get_nodes() const -> const auto & {return p_nodes;};
    inline auto get_rotations() const -> const auto & {return rotations;};
    inline auto get_center() const {return Base::world_coordinates(LocalCoordinates(0));};
    inline auto get_contains_local(const LocalCoordinates & point_i, const FLOATING_POINT_TYPE & eps) const -> bool {
        const auto & u = point_i[0];
        const auto & v = point_i[1];
        const auto & w = point_i[2];

        return IN_CLOSED_INTERVAL(this->p_nodes.row(0)[0] - eps, u, this->p_nodes.row(1)[0]) &&
                IN_CLOSED_INTERVAL(this->p_nodes.row(0)[1] - eps, v, this->p_nodes.row(1)[1]) &&
                IN_CLOSED_INTERVAL(this->p_nodes.row(0)[2] - eps, w, this->p_nodes.row(1)[2]);
    }


    friend struct Element<Beam>
    inline auto get_L(const LocalCoordinates & point_i) const -> Vector<NumberOfNodesAtCompileTime> {
        const auto  & u = point_i[0];
        const auto & v = point_i[1];
        const auto & w = point_i[2]; 


        return {
            /* à completer avec la matrice 3x3 des fonctions 
             de forme en fonction de u, v et w */
        };
    };

    inline auto get_dL(const LocalCoordinates & /*xi*/) const -> Matrix<NumberOfNodesAtCompileTime, CanonicalDimension> {
        return {
            /* à completer avec la matrice 3x3 des fonctions 
             de forme en fonction de u, v et w */
        };
    };

    inline auto get_gauss_nodes() const -> const auto & {
        static std::vector<GaussNode> gauss_nodes {
            GaussNode {LocalCoordinates(0), 2} // Node 0
        };
        return gauss_nodes;
    }

    
protected:
    Matrix<NumberOfNodesAtCompileTime, Dimension> p_nodes;
    Matrix<NumberOfNodesAtCompileTime, Dimension> rotations;

}

} // caribou::geometry