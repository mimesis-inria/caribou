cmake_minimum_required(VERSION 3.12)
project(Caribou)

# Manually define VERSION
set(Caribou_VERSION_MAJOR 1)
set(Caribou_VERSION_MINOR 0)
set(Caribou_VERSION_PATCH 0)
set(Caribou_VERSION ${Caribou_VERSION_MAJOR}.${Caribou_VERSION_MINOR}.${Caribou_VERSION_PATCH})

# Policies
cmake_policy(SET CMP0072 NEW)
cmake_policy(SET CMP0074 NEW)

## Default build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "RelWithDebInfo")
endif()

# Include CMake specific tools
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/CMake")
include(CheckCXXCompilerFlag)
include(CMakePackageConfigHelpers)
include(CMakeDependentOption)
include(caribouTools)
set(CMAKE_REQUIRED_QUIET 1)

CHECK_CXX_COMPILER_FLAG("-march=native" COMPILER_SUPPORTS_MARCH_NATIVE)

# OPTIONS
option(CARIBOU_USE_DOUBLE "Use double precision for floating point values." ON)
option(CARIBOU_BUILD_TESTS "Compile the automatic tests for Caribou." OFF)
CMAKE_DEPENDENT_OPTION(CARIBOU_OPTIMIZE_FOR_NATIVE "Build with -march=native to allow msse support" ON "COMPILER_SUPPORTS_MARCH_NATIVE" OFF)

# Tests
if (CARIBOU_BUILD_TESTS)
    # SOFA is exporting its own gtest, let's try to find it first
    find_package(GTest QUIET)
    if (NOT GTest_FOUND)
        # This find_package needs to be executed at MOST once, else it will fails on the second call
        find_package(GMock QUIET REQUIRED)
    endif ()
endif()

# SOFA option
find_package(SOFA COMPONENTS SofaFramework QUIET)
CMAKE_DEPENDENT_OPTION(CARIBOU_WITH_SOFA "Compile the SofaCaribou library" ON "SOFA_FOUND" OFF)
if (CARIBOU_WITH_SOFA)
    message(STATUS "Caribou with SOFA support\n\tSOFA root: ${SOFA_ROOT}")
endif()

# Python 3 option
set(PYBIND11_PYTHON_VERSION 3)
find_package(pybind11 CONFIG QUIET)
find_package(Python3 QUIET)
CMAKE_DEPENDENT_OPTION(CARIBOU_WITH_PYTHON_3 "Compile the plugin with python 3 support" ON "Python3_FOUND;pybind11_FOUND" OFF)

# MKL option
set(BUILD_SHARED_LIBS ON)
set(BLA_VENDOR Intel10_64lp)
find_package(BLAS QUIET)
find_package(LAPACKE QUIET)
find_package(MKL QUIET)
CMAKE_DEPENDENT_OPTION(CARIBOU_WITH_MKL "Compile the plugin with MKL support (uses LAPACK(E) and BLAS)" ON "MKL_FOUND;LAPACKE_FOUND;BLAS_FOUND" OFF)
if (CARIBOU_WITH_MKL)
    set(CARIBOU_MKL_LIBRARIES "${MKL_LIBRARIES};${LAPACKE_LIBRARIES};${BLAS_LIBRARIES}")
    list(REMOVE_DUPLICATES CARIBOU_MKL_LIBRARIES)
    message(STATUS "Caribou with MKL support\n\tInclude directories: ${MKL_INCLUDE_DIRS}\n\tLibraries: ${CARIBOU_MKL_LIBRARIES}")
endif()

# OpenMP option
find_package(OpenMP QUIET)
CMAKE_DEPENDENT_OPTION(CARIBOU_WITH_OPENMP "Compile the plugin with OpenMP support" ON "OPENMP_FOUND" OFF)
if (CARIBOU_WITH_OPENMP)
    message(STATUS "Caribou with OpenMP support\n\tLibraries: ${OpenMP_CXX_LIBRARIES}")
endif()


## Change default install prefix
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX ${CMAKE_BINARY_DIR}/install CACHE PATH "Install path prefix, prepended onto install directories." FORCE)
endif()
message(STATUS "Installation prefix: ${CMAKE_INSTALL_PREFIX}")

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Compiler options
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED 17)
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
