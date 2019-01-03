cmake_minimum_required(VERSION 3.1)
project(Caribou)

## Default build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "RelWithDebInfo")
endif()


# OPTIONS

option(CARIBOU_BUILD_TESTS "Compile the automatic tests for Caribou." ON)
option(CARIBOU_WITH_SOFA "Compile the SofaCaribou library" ON)

find_package(PythonLibs 2 QUIET)
option(CARIBOU_WITH_PYTHON_2 "Compile the plugin with python 2 support" ${PYTHONLIBS_FOUND})

find_package(PythonLibs 3 QUIET)
option(CARIBOU_WITH_PYTHON_3 "Compile the plugin with python 3 support" ${PYTHONLIBS_FOUND})


## Change default install prefix
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX ${CMAKE_BINARY_DIR}/install CACHE PATH "Install path prefix, prepended onto install directories." FORCE)
endif()
message(STATUS "Install prefix: ${CMAKE_INSTALL_PREFIX}")

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
endif()
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic")

add_subdirectory(src/Caribou)

if (CARIBOU_WITH_SOFA)
    add_subdirectory(src/SofaCaribou)
endif()