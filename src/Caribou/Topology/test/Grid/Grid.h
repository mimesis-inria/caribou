#ifndef CARIBOU_TOPOLOGY_TEST_GRID_H
#define CARIBOU_TOPOLOGY_TEST_GRID_H
#include <Eigen/Core>
#include <list>


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

#include "Grid_1D.h"
#include "Grid_2D.h"
#include "Grid_3D.h"
#endif //CARIBOU_TOPOLOGY_TEST_GRID_H
