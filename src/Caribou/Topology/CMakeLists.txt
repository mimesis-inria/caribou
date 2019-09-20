cmake_minimum_required(VERSION 3.1)
project(Caribou.Topology)

set(HEADER_FILES
    Grid/Grid.h
    Grid/Internal/BaseGrid.h
    Grid/Internal/BaseMultidimensionalGrid.h
    Grid/Internal/BaseUnidimensionalGrid.h
    HashGrid.h)

add_library(${PROJECT_NAME} INTERFACE)

target_include_directories(${PROJECT_NAME} INTERFACE "$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src/>")
target_include_directories(${PROJECT_NAME} INTERFACE "$<INSTALL_INTERFACE:include>")

target_link_libraries(${PROJECT_NAME} INTERFACE Caribou.Geometry)

if (CARIBOU_BUILD_TESTS)
    add_subdirectory(test)
endif()

if (CARIBOU_WITH_PYTHON_2 OR CARIBOU_WITH_PYTHON_3)
    add_subdirectory(Python)
endif()