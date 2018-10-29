cmake_minimum_required(VERSION 3.1)
project(Caribou)

add_subdirectory(Components)
add_subdirectory(Algebra)
add_subdirectory(Geometry)
add_subdirectory(Topology)
add_subdirectory(Plugin)

if (CARIBOU_WITH_PYTHON_2 OR CARIBOU_WITH_PYTHON_3)
    add_subdirectory(Python)
endif()