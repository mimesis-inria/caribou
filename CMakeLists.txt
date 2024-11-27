cmake_minimum_required(VERSION 3.12)
project(Caribou VERSION 21.06.00)

# Policies
cmake_policy(SET CMP0072 NEW)
cmake_policy(SET CMP0074 NEW)
cmake_policy(SET CMP0076 NEW)
cmake_policy(SET CMP0094 NEW)

set (CARIBOU_VERSION ${PROJECT_VERSION})

## Default build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()

## Change default install prefix
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX ${CMAKE_BINARY_DIR}/install CACHE PATH "Install path prefix, prepended onto install directories." FORCE)
endif()

# Detect if Caribou is a subproject of SOFA
if ("${CMAKE_PROJECT_NAME}" STREQUAL "Sofa")
    set(CARIBOU_COMPILED_AS_SOFA_SUBPROJECT 1)
    message(STATUS "${PROJECT_NAME} is compiled as a subproject to SOFA.")
endif ()

# If the plugin is compiled within SOFA, add an option to deactivate it, in which case the remaining
# of this file won't run and therefore will not generate any error if the requirements (python, pybind11)
# are not found.
if (CARIBOU_COMPILED_AS_SOFA_SUBPROJECT)
    option(PLUGIN_SOFACARIBOU "Build the caribou plugin." ON)
    if (NOT PLUGIN_SOFACARIBOU)
        return()
    else()
        message("Adding plugin caribou")
    endif()
endif()

# Include CMake specific tools
list(INSERT CMAKE_MODULE_PATH 0 "${CMAKE_CURRENT_LIST_DIR}/cmake")
list(INSERT CMAKE_PREFIX_PATH 0 "${CMAKE_BINARY_DIR}/cmake")
include(CheckCXXCompilerFlag)
include(CMakePackageConfigHelpers)
include(CMakeDependentOption)
include(caribouTools)
set(CMAKE_REQUIRED_QUIET 1)

CHECK_CXX_COMPILER_FLAG("-march=native" COMPILER_SUPPORTS_MARCH_NATIVE)

# OPTIONS
option(CARIBOU_USE_FLOAT "Use float precision for floating point values." OFF)
CMAKE_DEPENDENT_OPTION(CARIBOU_OPTIMIZE_FOR_NATIVE "Build with -march=native to allow msse support" ON "COMPILER_SUPPORTS_MARCH_NATIVE" OFF)

# SOFA option
find_package(Threads QUIET)
find_package(SOFA COMPONENTS Sofa.Config QUIET)
CMAKE_DEPENDENT_OPTION(CARIBOU_WITH_SOFA "Compile the SofaCaribou library" ON "SOFA_FOUND" OFF)
if (CARIBOU_WITH_SOFA)
    list(APPEND CMAKE_PREFIX_PATH ${SOFA_MODULE_DIR})
    list(REMOVE_DUPLICATES CMAKE_PREFIX_PATH)
    find_package(SofaPython3 QUIET) # We need to do this first to setup the exact python version required
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
        message(STATUS "Caribou with unittests")
    else()
        add_library(gtest SHARED IMPORTED)
        set_target_properties(gtest PROPERTIES
            INTERFACE_INCLUDE_DIRECTORIES "${GTEST_INCLUDE_DIRS}"
            IMPORTED_LOCATION "${GTEST_BOTH_LIBRARIES}"
        )
        message(STATUS "Caribou with unittests\n\tLibraries: ${GTEST_BOTH_LIBRARIES}\n\tInclude directories: ${GTEST_INCLUDE_DIRS}")
    endif ()
elseif(NOT GTest_FOUND AND NOT GMock_FOUND)
    message(STATUS "Caribou without unittests (missing GTest)")
endif()

# Benchmarks
OPTION(CARIBOU_BUILD_BENCHMARKS "Compile the benchmarks for Caribou." ON)
if (CARIBOU_BUILD_BENCHMARKS)
    message(STATUS "Caribou with benchmarks")
else()
    message(STATUS "Caribou without benchmarks")
endif()

# Python 3 option
if (CARIBOU_WITH_SOFA)
    sofa_find_package(SofaPython3 QUIET)
endif()
set(CMAKE_CXX_STANDARD 17)
set(PYBIND11_FINDPYTHON ON)
set(PYBIND11_PYTHON_VERSION 3)
set(Python_FIND_FRAMEWORK LAST)
find_package(SofaPython3 QUIET)
if (SofaPython3_FOUND)
    set(PYTHON_EXECUTABLE ${PYTHON_EXECUTABLE})
    set(PYBIND11_PYTHON_VERSION ${PYTHON_VERSION})
endif()
find_package(Python ${PYBIND11_PYTHON_VERSION} COMPONENTS Interpreter Development QUIET)
find_package(pybind11 CONFIG QUIET)
CMAKE_DEPENDENT_OPTION(CARIBOU_WITH_PYTHON_3 "Compile the plugin with python 3 support" ON "Python_FOUND;pybind11_FOUND" OFF)
if (CARIBOU_WITH_PYTHON_3)
    message(STATUS "Caribou with Python support:
    Version: ${Python_VERSION}
    Executable: ${Python_EXECUTABLE}
    Headers: ${Python_INCLUDE_DIRS}
    Libraries: ${Python_LIBRARIES}
    pybind11 version: ${pybind11_VERSION}
    pybind11 config: ${pybind11_CONFIG}"
)
endif()

# SOFA Python 3
if (CARIBOU_WITH_SOFA)
    get_filename_component(SOFA_ROOT_DIR "${SofaFramework_DIR}/../../.." ABSOLUTE)
    message(STATUS
            "Caribou with SOFA support\n"
            "    SOFA version: ${SofaFramework_VERSION}\n"
            "    SOFA location: ${SOFA_ROOT_DIR}"
            )
    if (CARIBOU_WITH_PYTHON_3)
        if (SofaPython3_FOUND)
            get_filename_component(SP3_ROOT_DIR "${SofaPython3_DIR}/../../.." ABSOLUTE)
            message(STATUS
                    "Caribou with SofaPython3 support\n"
                    "    SP3 version: ${SofaPython3_VERSION}\n"
                    "    SP3 location: ${SP3_ROOT_DIR}"
            )
        else()
            message(STATUS "Caribou without SofaPython (not found)")
        endif()
    endif ()
endif()

# OpenMP option
set(CMAKE_FIND_LIBRARY_SUFFIXES_BEFORE ${CMAKE_FIND_LIBRARY_SUFFIXES})
SET(CMAKE_FIND_LIBRARY_SUFFIXES .so)
find_package(OpenMP QUIET)
set(CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES_BEFORE})
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
set(CARIBOU_VTK_MODULES CommonCore CommonDataModel IOLegacy)
find_package(VTK QUIET)
if (VTK_FOUND AND VTK_VERSION VERSION_LESS "9.0")
    set(CARIBOU_VTK_MODULES vtkCommonCore vtkCommonDataModel vtkIOLegacy)
endif()
find_package(VTK COMPONENTS ${CARIBOU_VTK_MODULES} QUIET)
CMAKE_DEPENDENT_OPTION(CARIBOU_WITH_VTK "Compile the plugin with VTK support." ON "VTK_FOUND" OFF)
if (CARIBOU_WITH_VTK)
    foreach(VTK_MODULE ${CARIBOU_VTK_MODULES})
        list(APPEND _VTK_INCLUDE_DIRECTORIES "${${VTK_MODULE}_INCLUDE_DIRS}")
        list(APPEND _VTK_LIBRARY_DIRECTORIES "${${VTK_MODULE}_RUNTIME_LIBRARY_DIRS}")
        list(APPEND _VTK_LIBRARIES "${${VTK_MODULE}_LIBRARIES}")
    endforeach()
    list(REMOVE_DUPLICATES _VTK_INCLUDE_DIRECTORIES)
    list(REMOVE_DUPLICATES _VTK_LIBRARY_DIRECTORIES)
    list(REMOVE_DUPLICATES _VTK_LIBRARIES)
    set(VTK_MSG "Caribou with VTK support")
    if (VTK_PREFIX_PATH)
        set(VTK_MSG "${VTK_MSG}\n\tPrefix: ${VTK_PREFIX_PATH}")
    endif()
    if (_VTK_INCLUDE_DIRECTORIES)
        set(VTK_MSG "${VTK_MSG}\n\tInclude directories: ${_VTK_INCLUDE_DIRECTORIES}")
    endif ()
    if (_VTK_LIBRARY_DIRECTORIES)
        set(VTK_MSG "${VTK_MSG}\n\tLibrary directories: ${_VTK_LIBRARY_DIRECTORIES}")
    endif ()
    if (_VTK_LIBRARIES)
        set(VTK_MSG "${VTK_MSG}\n\tLibraries directories: ${_VTK_LIBRARIES}")
    endif ()
    set(VTK_MSG "${VTK_MSG}\n\tVersion: ${VTK_VERSION}")
    message(STATUS "${VTK_MSG}")
endif()

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Compiler options
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED 17)
if(NOT WIN32)
    if ("${CMAKE_CXX_COMPILER} ${CMAKE_CXX_COMPILER_ARG1}" MATCHES ".*[cC][lL][aA][nN][gG]")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fcolor-diagnostics")
    else()
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fdiagnostics-color")
    endif()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic")
endif()

add_subdirectory(src/Caribou)

if (CARIBOU_WITH_SOFA)
    add_subdirectory(src/SofaCaribou)
endif()

if (CARIBOU_BUILD_TESTS)
    add_subdirectory(unittest)
endif()

if (CARIBOU_BUILD_BENCHMARKS)
    add_subdirectory(Benchmark)
endif()

install(
    FILES
        "${CMAKE_CURRENT_SOURCE_DIR}/cmake/caribouTools.cmake"
    DESTINATION lib/cmake/Caribou
    COMPONENT headers
)

install(
        FILES
        "${CMAKE_CURRENT_SOURCE_DIR}/cmake/caribouTools.cmake"
        "${CMAKE_CURRENT_SOURCE_DIR}/cmake/FindSOFA.cmake"
        "${CMAKE_CURRENT_SOURCE_DIR}/cmake/FindMKL.cmake"
        DESTINATION lib/cmake/SofaCaribou
        COMPONENT headers
)
