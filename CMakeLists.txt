cmake_minimum_required(VERSION 3.5)
project(Caribou)

# Manually define VERSION
set(Caribou_VERSION_MAJOR 1)
set(Caribou_VERSION_MINOR 0)
set(Caribou_VERSION_PATCH 0)
set(Caribou_VERSION ${Caribou_VERSION_MAJOR}.${Caribou_VERSION_MINOR}.${Caribou_VERSION_PATCH})

## Default build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "RelWithDebInfo")
endif()

# CMAKE TOOLS
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/CMake")
include(caribouTools)

# OPTIONS
option(CARIBOU_USE_DOUBLE "Use double precision for floating point values." ON)
option(CARIBOU_BUILD_TESTS "Compile the automatic tests for Caribou." OFF)
option(CARIBOU_WITH_SOFA "Compile the SofaCaribou library" ON)

find_package(pybind11 CONFIG QUIET)

# Python 2 option
find_package(PythonLibs 2 QUIET)
if (${PYTHONLIBS_FOUND} AND ${pybind11_FOUND})
    option(CARIBOU_WITH_PYTHON_2 "Compile the plugin with python 2 support" ON)
else()
    option(CARIBOU_WITH_PYTHON_2 "Compile the plugin with python 2 support" OFF)
endif()

# Python 3 option
find_package(Python3 QUIET)
if (${Python3_FOUND} AND ${pybind11_FOUND})
    option(CARIBOU_WITH_PYTHON_3 "Compile the plugin with python 3 support" ON)
else()
    option(CARIBOU_WITH_PYTHON_3 "Compile the plugin with python 3 support" OFF)
endif()

# OpenMP option
find_package(OpenMP)
if (OPENMP_FOUND)
    OPTION(CARIBOU_WITH_OPENMP "Use OpenMP support" ON)
else()
    OPTION(CARIBOU_WITH_OPENMP "Use OpenMP support" OFF)
endif()
if (CARIBOU_WITH_OPENMP)
    find_package(OpenMP REQUIRED)
endif()


## Change default install prefix
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX ${CMAKE_BINARY_DIR}/install CACHE PATH "Install path prefix, prepended onto install directories." FORCE)
endif()
message(STATUS "Install prefix: ${CMAKE_INSTALL_PREFIX}")

include(CMakePackageConfigHelpers)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

message(STATUS "Build prefix: ${CMAKE_BINARY_DIR}")
message(STATUS "Archive prefix: ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}")
message(STATUS "Library prefix: ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}")
message(STATUS "Runtime prefix: ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")

# Compiler options
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")
if ("${CMAKE_CXX_COMPILER} ${CMAKE_CXX_COMPILER_ARG1}" MATCHES ".*[cC][lL][aA][nN][gG]")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fcolor-diagnostics")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fdiagnostics-color")
endif()
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic")

add_subdirectory(src/Caribou)

if (CARIBOU_WITH_SOFA)
    add_subdirectory(src/SofaCaribou)
endif()
