#ifndef CARIBOU_GEOMETRY_NODE_H
#define CARIBOU_GEOMETRY_NODE_H

#include <cassert>
#include <numeric>
#include <cstddef>
#include <array>
#include <initializer_list>
#include <algorithm>

#include <Caribou/config.h>
#include <Caribou/Geometry/Internal/BaseNode.h>
#include <Caribou/Algebra/Vector.h>

namespace caribou
{
namespace geometry
{

template<size_t Dim, REQUIRES(Dim > 0 and Dim < 4)>
class Node : public internal::BaseNode<Dim, Node<Dim>>
{
};

/**
 * A 1D node in space.
 */
 template <>
class Node<1>  : public internal::BaseNode<1, Node<1>>
{
    using Self = Node<1>;
    using Base = BaseNode<1, Node<1>>;

public:
    static constexpr size_t Dimension = 1;
    using ValueType = typename Base::ValueType;

    using Base::Base; // Import constructors

    /**
     * Forwarding constructor
     * This constructor can be used to construct a Node from any class that does not inherits Node
     * via specialisation of the OtherNodeType template.
     */
    template <
            typename AnyType,
            REQUIRES(not std::is_base_of_v<algebra::internal::CaribouMatrix, AnyType>)
    >
    Node(AnyType && anything) : Base (anything) {

    }

    /** Assignment operator **/
    template<typename OtherValueType>
    Node &
    operator = (const caribou::algebra::Vector<Dimension, OtherValueType> & v)
    {
        for (size_t i = 0; i < Dimension; ++i)
            (*this)[i] = v[i];
        return (*this);
    }

    inline const ValueType & x () const { return this->at(0); }
    inline ValueType & x () { return this->at(0); }
};

/**
 * A 2D node in space.
 */
template<>
class Node<2>  : public internal::BaseNode<2, Node<2>>
{
    using Self = Node<2>;
    using Base = internal::BaseNode<2, Node<2>>;

public:
    static constexpr size_t Dimension = 2;
    using ValueType = typename Base::ValueType;

    using Base::Base; // Import constructors

    /**
     * Forwarding constructor
     * This constructor can be used to construct a Node from any class that does not inherits Node
     * via specialisation of the OtherNodeType template.
     */
    template <
            typename AnyType,
            REQUIRES(not std::is_base_of_v<algebra::internal::CaribouMatrix, AnyType>)
    >
    Node(AnyType && anything) : Base (anything) {

    }

    /** Assignment operator **/
    template<typename OtherValueType>
    Node &
    operator = (const caribou::algebra::Vector<Dimension, OtherValueType> & v)
    {
        for (size_t i = 0; i < Dimension; ++i)
            (*this)[i] = v[i];
        return (*this);
    }

    inline const ValueType & x () const { return this->at(0); }
    inline ValueType & x () { return this->at(0); }

    inline const ValueType & y () const { return this->at(1); }
    inline ValueType & y () { return this->at(1); }
};

/**
 * A 3D node in space.
 */
template<>
class Node<3>  : public internal::BaseNode<3, Node<3>>
{
    using Self = Node<3>;
    using Base = internal::BaseNode<3, Node<3>>;

public:
    static constexpr size_t Dimension = 3;
    using ValueType = typename Base::ValueType;

    using Base::Base;

    /**
     * Forwarding constructor
     * This constructor can be used to construct a Node from any class that does not inherits Node
     * via specialisation of the OtherNodeType template.
     */
    template <
            typename AnyType,
            REQUIRES(not std::is_base_of_v<algebra::internal::CaribouMatrix, AnyType>)
    >
    Node(AnyType && anything) : Base (anything) {

    }

    /** Assignment operator **/
    template<typename OtherValueType>
    Node &
    operator = (const caribou::algebra::Vector<Dimension, OtherValueType> & v)
    {
        for (size_t i = 0; i < Dimension; ++i)
            (*this)[i] = v[i];
        return (*this);
    }

    inline const ValueType & x () const { return this->at(0); }
    inline ValueType & x () { return this->at(0); }

    inline const ValueType & y () const { return this->at(1); }
    inline ValueType & y () { return this->at(1); }

    inline const ValueType & z () const { return this->at(2); }
    inline ValueType & z () { return this->at(2); }
};

// Deduction guides

/**
 * Node(x) or Node(x,y) or Node(x,y,z)
 * where x,y and z are of the same type which is integral or floating point.
 */
template <
        typename T, typename... Ts,
        REQUIRES(std::is_integral_v<T> or std::is_floating_point_v<T>),
        REQUIRES(std::conjunction_v<std::is_same<T, Ts>...>)
>
Node(T, Ts...) -> Node <sizeof...(Ts)+1>;

/**
 * Node(n)
 * where n is another Node type (type that derives BaseNode).
 */
template <
        typename T,
        REQUIRES(std::is_base_of_v<internal::BaseNode, T>)
>
Node(T) -> Node <T::Dimension>;

/**
 * Node(c)
 * where c are the coordinates. The type of c must be an array type (c-array, eg. T[Dim]) of rank 1 (eg. not T[][]) and
 * the value of the coordinates must be convertible to floating point type.
 */
template <
        typename ValueType,
        size_t Dim,
        /* ValueType must be convertible to floating points */
        REQUIRES(std::is_convertible_v<ValueType, FLOATING_POINT_TYPE>)
>
Node(const ValueType (&) [Dim]) -> Node <Dim>;

/**
* Node(c)
* where c is a std::array that contains the coordinates.
*/
template <
        typename T,
        std::size_t Dim,
        /* CoordinateType must contain elements convertible to floating points */
        REQUIRES(std::is_convertible_v<T, FLOATING_POINT_TYPE>)
>
Node(const std::array<T, Dim> &) -> Node <Dim>;

} // namespace geometry

} // namespace caribou
#endif //CARIBOU_GEOMETRY_NODE_H
