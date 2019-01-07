#ifndef CARIBOU_GEOMETRY_INTERNAL_BASENODE_H
#define CARIBOU_GEOMETRY_INTERNAL_BASENODE_H

#include <Caribou/config.h>
#include <Caribou/Algebra/Vector.h>

namespace caribou {
namespace geometry {
namespace internal {

/**
 * A node in space (independent of the space dimension).
 *
 * ** Do not use this class directly. Use instead caribou::algebra::Matrix. **
 * @tparam Dim Dimension of the current space.
 * @tparam NodeType_ The derived node class
 */
template<
        size_t Dim,
        typename NodeType_
>
class BaseNode : public caribou::algebra::Vector<Dim, FLOATING_POINT_TYPE>
{
    static_assert(Dim > 0 and Dim < 4, "Only nodes of dimension 1, 2 or 3 are permitted.");

public:
    static constexpr size_t Dimension = Dim;
    using VectorType = caribou::algebra::Vector<Dimension, FLOATING_POINT_TYPE>;
    using ValueType = typename VectorType::ValueType;
    using NodeType = NodeType_;

    using VectorType::VectorType;

    /** Default constructor **/
    constexpr BaseNode() : VectorType() {}

    /** Default constructor with zero initialization **/
    constexpr BaseNode(bool initialize_to_zero) : VectorType(initialize_to_zero) {}

    /** Copy constructor by vector type **/
    constexpr BaseNode(const VectorType & v) : VectorType(v) {}

    /** Copy constructor from another Node **/
    template<
            typename OtherNodeType,
            REQUIRES(std::is_base_of_v<algebra::internal::CaribouMatrix, OtherNodeType>)
    >
    constexpr BaseNode(const OtherNodeType & p) : VectorType (p) {
        static_assert(OtherNodeType::Dimension == Dimension, "Cannot construct a Node from another Node of different dimension.");
    }

    /**
     * Forwarding constructor
     * This constructor can be used to construct a Node from any class that does not inherits Node
     * via specialisation of the OtherNodeType template.
     */
    template <
            typename AnyType,
            REQUIRES(not std::is_base_of_v<algebra::internal::CaribouMatrix, AnyType>)
    >
    BaseNode(AnyType && anything) : VectorType (anything) {

    }

    /**
     * Return a copy of this node translated by the vector t
     */
    inline NodeType
    translated(const VectorType & t) const
    {
        return NodeType(*this).translate(t);
    }


    /**
     * Translate this node by the vector t
     */
    inline NodeType &
    translate(const VectorType & t)
    {
        for (size_t i = 0; i < Dimension; ++i) {
            (*this)[i] += t[i];
        }

        return static_cast<NodeType&>(*this);
    }
};

} // namespace internal
} // namespace geometry
} // namespace caribou
#endif //CARIBOU_GEOMETRY_INTERNAL_BASENODE_H
