#ifndef CARIBOU_TOPOLOGY_TEST_GRID_H
#define CARIBOU_TOPOLOGY_TEST_GRID_H
#include <Eigen/Core>
#include <list>
#include <Caribou/Geometry/Traits.h>

template<template<typename, typename...> typename H, typename T, typename... Ts>
bool list_are_equals (const std::list<T> & l1, const H<T, Ts...> & l2)
{
    if (l1.size() != l2.size())
        return false;

    typename std::list<T>::const_iterator it1;
    typename H<T, Ts...>::const_iterator it2;

    for (it1 = l1.begin(), it2 = l2.begin(); it1 != l1.end(); ++it1, ++it2)
        if (*it1 != *it2)
            return false;
    return true;
}

template <typename Node>
void EXPECT_NODE_EQ(const Node & n1, const Node & n2) {
    EXPECT_FLOAT_EQ((n2-n1).norm(), 0.);
}

template <typename Element>
void EXPECT_ELEMENTS_EQ(const Element & e1, const Element & e2) {
    for (std::size_t i = 0; i < Element::NumberOfNodes; ++i) {
        EXPECT_NODE_EQ(e1.node(i), e2.node(i));
    }
}

template <typename Element1, typename Element2>
void EXPECT_ELEMENTS_EQ(const Element1 & e1, const Element2 & e2) {
    static_assert((int)caribou::traits<Element1>::NumberOfNodes == (int)caribou::traits<Element2>::NumberOfNodes);
    for (std::size_t i = 0; i < Element1::NumberOfNodes; ++i) {
        EXPECT_NODE_EQ(e1.node(i), e2.node(i));
    }
}

template <typename Element, typename... Elements>
void EXPECT_ELEMENTS_EQ(const Element & e1, const Element & e2, Elements... elements) {
    EXPECT_ELEMENTS_EQ(e1, e2);
    EXPECT_ELEMENTS_EQ(e2, std::forward<Elements>(elements)...);
}

#include "Grid_1D.h"
#include "Grid_2D.h"
#include "Grid_3D.h"
#endif //CARIBOU_TOPOLOGY_TEST_GRID_H
