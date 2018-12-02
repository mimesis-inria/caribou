cmake_minimum_required(VERSION 3.1)

option(CARIBOU_USE_DOUBLE "Use double precision for floating point values." OFF)

add_subdirectory(Algebra)
add_subdirectory(Geometry)
add_subdirectory(Mechanics)
add_subdirectory(Topology)

if (CARIBOU_WITH_PYTHON_2 OR CARIBOU_WITH_PYTHON_3)
    add_subdirectory(Python)
endif()