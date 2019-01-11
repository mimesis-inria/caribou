#ifndef SOFACARIBOU_TRAITS_H
#define SOFACARIBOU_TRAITS_H

#include <sofa/defaulttype/Vec.h>
#include <Caribou/Geometry/Node.h>

namespace caribou::geometry {

//
// sofa Vec -> caribou Node
//

/**
* Node(v)
* where v is a sofa vector that contains the coordinates.
*/
template <
        typename T,
        std::size_t Dim
>
Node(const sofa::defaulttype::Vec<Dim, T> &) -> Node <Dim>;
template <>
Node<1>::Node(const sofa::defaulttype::Vec<1, float> & v) : Node<1>(v[0]) {}
template <>
Node<1>::Node(const sofa::defaulttype::Vec<1, double> & v) : Node<1>(v[0]) {}
template <>
Node<2>::Node(const sofa::defaulttype::Vec<2, float> & v) : Node<2>(v[0], v[1]) {}
template <>
Node<2>::Node(const sofa::defaulttype::Vec<2, double> & v) : Node<2>(v[0], v[1]) {}
template <>
Node<3>::Node(const sofa::defaulttype::Vec<3, float> & v) : Node<3>(v[0], v[1], v[2]) {}
template <>
Node<3>::Node(const sofa::defaulttype::Vec<3, double> & v) : Node<3>(v[0], v[1], v[2]) {}
}

#endif //SOFACARIBOU_TRAITS_H
