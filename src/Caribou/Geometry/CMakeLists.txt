cmake_minimum_required(VERSION 3.1)
project(Caribou.Geometry)

set(HEADER_FILES
    Hexahedron.h
    Node.h
    Quad.h
    Segment.h
    Triangle.h
    Internal/BaseNode.h
    Internal/BaseQuad.h
    Internal/BaseTriangle.h
    Interpolation/CanonicalElement.h
    Interpolation/Hexahedron.h
    Interpolation/Quad.h
    Interpolation/Triangle.h)

add_library(${PROJECT_NAME} INTERFACE)

target_include_directories(${PROJECT_NAME} INTERFACE "$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src/>")
target_include_directories(${PROJECT_NAME} INTERFACE "$<INSTALL_INTERFACE:include>")

if (CARIBOU_BUILD_TESTS)
    add_subdirectory(test)
endif()

if (CARIBOU_WITH_PYTHON_2 OR CARIBOU_WITH_PYTHON_3)
    add_subdirectory(Python)
endif()