#ifndef SOFACARIBOU_TRAITS_H
#define SOFACARIBOU_TRAITS_H

#include <sofa/defaulttype/Vec.h>
#include <Caribou/Geometry/Node.h>

namespace caribou::geometry {

/**
* Node(v)
* where v is a sofa vector that contains the coordinates.
*/
template <
        typename T,
        std::size_t Dim
>
Node(const sofa::defaulttype::Vec<Dim, T> &) -> Node <Dim>;


template<>
Node<3>::Node(const sofa::defaulttype::Vec<3, double> & v) : Node {v[0], v[1], v[2]} {}

//template <
//        typename T,
//        std::size_t Dim
//>
//struct internal::NodeConverter<const sofa::defaulttype::Vec<Dim, T> &> {
//    static inline Node<Dim> convert(const sofa::defaulttype::Vec<Dim, T> & v) {
//        Node<Dim> n;
//        for (size_t i = 0; i < Dim; ++i)
//            n[i] = v[i];
//        return n;
//    };
//};

}

#endif //SOFACARIBOU_TRAITS_H
