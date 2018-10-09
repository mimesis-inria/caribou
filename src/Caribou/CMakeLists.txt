cmake_minimum_required(VERSION 3.1)
project(Caribou)

add_subdirectory(Geometry)

if (CARIBOU_WITH_PYTHON_2 OR CARIBOU_WITH_PYTHON_3)
    add_subdirectory(Python)
endif()