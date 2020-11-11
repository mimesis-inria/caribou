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
CMAKE_DEPENDENT_OPTION(CARIBOU_OPTIMIZE_FOR_NATIVE "Build with -march=native to allow msse support" ON "COMPILER_SUPPORTS_MARCH_NATIVE" OFF)

# SOFA option
find_package(SOFA COMPONENTS SofaFramework QUIET)
CMAKE_DEPENDENT_OPTION(CARIBOU_WITH_SOFA "Compile the SofaCaribou library" ON "SOFA_FOUND" OFF)
if (CARIBOU_WITH_SOFA)
    list(APPEND CMAKE_PREFIX_PATH ${SOFA_MODULE_DIR})
    list(REMOVE_DUPLICATES CMAKE_PREFIX_PATH)
    message(STATUS "${SOFA_MODULE_DIR}")
endif()

# Tests
# SOFA is exporting its own gtest, let's try to find it first
find_package(GTest CONFIG QUIET)
if (NOT GTest_FOUND)
    # This find_package needs to be executed at MOST once, else it will fails on the second call
    find_package(GMock QUIET)
endif ()
CMAKE_DEPENDENT_OPTION(CARIBOU_BUILD_TESTS "Compile the automatic tests for Caribou." ON "GTest_FOUND OR GMock_FOUND" OFF)
if (CARIBOU_BUILD_TESTS)
    if (TARGET gtest)
        get_target_property(GTEST_INCLUDE_DIRS gtest INTERFACE_INCLUDE_DIRECTORIES)
        get_target_property(GTEST_BOTH_LIBRARIES gtest LOCATION)
    else()
        add_library(gtest SHARED IMPORTED)
        set_target_properties(gtest PROPERTIES
            INTERFACE_INCLUDE_DIRECTORIES "${GTEST_INCLUDE_DIRS}"
            IMPORTED_LOCATION "${GTEST_BOTH_LIBRARIES}"
        )
    endif ()
    message(STATUS "Caribou with unittests\n\tLibraries: ${GTEST_BOTH_LIBRARIES}\n\tInclude directories: ${GTEST_INCLUDE_DIRS}")
endif()

# Python 3 option
set(PYBIND11_PYTHON_VERSION 3)
find_package(pybind11 CONFIG QUIET)
find_package(Python3 QUIET)
CMAKE_DEPENDENT_OPTION(CARIBOU_WITH_PYTHON_3 "Compile the plugin with python 3 support" ON "Python3_FOUND;pybind11_FOUND" OFF)
if (CARIBOU_WITH_PYTHON_3)
    message(STATUS "Caribou with Python support\n\tPython version: ${PYTHON_VERSION_STRING}\n\tPython executabe: ${PYTHON_EXECUTABLE}\n\tpybind11 version: ${pybind11_VERSION}\n\tpybind11 config: ${pybind11_CONFIG}")
endif()

# SOFA Python 3
if (CARIBOU_WITH_SOFA)
    find_package(SofaPython3 QUIET)
    if (CARIBOU_WITH_PYTHON_3 AND SofaPython3_FOUND)
        message(STATUS "Caribou with SOFA support\n\tSOFA: ${SOFA_ROOT}\n\tSofaPython3: ${SofaPython3_CONFIG}")
    else()
        message(STATUS "Caribou with SOFA support\n\tSOFA root: ${SOFA_ROOT}")
    endif ()
endif()

# OpenMP option
find_package(OpenMP QUIET)
CMAKE_DEPENDENT_OPTION(CARIBOU_WITH_OPENMP "Compile the plugin with OpenMP support" ON "OPENMP_FOUND" OFF)
if (CARIBOU_WITH_OPENMP)
    message(STATUS "Caribou with OpenMP support\n\tLibraries: ${OpenMP_CXX_LIBRARIES}")
endif()

# MKL option
set(MKL_STATIC ON)
if (NOT CARIBOU_WITH_OPENMP)
    set(MKL_THREADING_VENDOR SEQUENTIAL)
endif()
find_package(MKL QUIET)
CMAKE_DEPENDENT_OPTION(CARIBOU_WITH_MKL "Compile the plugin with MKL support." ON "MKL_FOUND" OFF)
if (CARIBOU_WITH_MKL)
    message(STATUS "Caribou with MKL support\n\tInclude directories: ${MKL_INCLUDE_DIRS}\n\tLibraries: ${MKL_LIBRARIES}\n\tFlags: ${MKL_LINKER_FLAGS}")
endif()

# VTK option
set(CARIBOU_VTK_MODULES CommonDataModel IOLegacy)
find_package(VTK QUIET)
if (VTK_FOUND AND VTK_VERSION VERSION_LESS "9.0")
    set(CARIBOU_VTK_MODULES vtkCommonDataModel vtkIOLegacy)
endif()
find_package(VTK COMPONENTS ${CARIBOU_VTK_MODULES} QUIET)
CMAKE_DEPENDENT_OPTION(CARIBOU_WITH_VTK "Compile the plugin with VTK support." ON "VTK_FOUND" OFF)
if (CARIBOU_WITH_VTK)
    message(STATUS "Caribou with VTK support\n\tVersion: ${VTK_VERSION}")
endif()

## Change default install prefix
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX ${CMAKE_BINARY_DIR}/install CACHE PATH "Install path prefix, prepended onto install directories." FORCE)
endif()

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
