cmake_minimum_required(VERSION 3.1)
project(caribou)

## Default build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "RelWithDebInfo")
endif()

# Try to find Sofa's cmake prefix
# First, make sure it isn't already been manually set by the user's cmake variables
find_path (SOFA_PREFIX SofaFramework/SofaFrameworkConfig.cmake)
if (SOFA_PREFIX)
    # The user already set the correct sofa cmake prefix
    message(STATUS "Sofa prefix already specified at path ${SOFA_PREFIX}")
else()
    # If it isn't found, let's try to find it from enviromnemnt variables such as SOFA_ROOT
    find_path (
            SOFA_PREFIX
            SofaFramework/SofaFrameworkConfig.cmake
            PATHS
            $ENV{SOFA_ROOT}/install/lib/cmake
            $ENV{SOFA_ROOT}/build/install/lib/cmake
            $ENV{SOFA_ROOT}/build/master/install/lib/cmake
            $ENV{SOFA_ROOT}/lib/cmake
            NO_DEFAULT_PATH
    )

    # If we found it, add it to cmake prefix path so that further find_package correctly locate sofa
    if (SOFA_PREFIX)
        message(STATUS "Sofa found at ${SOFA_PREFIX}")
        list(APPEND CMAKE_PREFIX_PATH "${SOFA_PREFIX}")
    endif()
endif()

## Change default install prefix
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX ${CMAKE_BINARY_DIR}/install CACHE PATH "Install path prefix, prepended onto install directories." FORCE)
endif()
message(STATUS "Install prefix: ${CMAKE_INSTALL_PREFIX}")

# Compiler options
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic")

add_subdirectory(plugin)
add_subdirectory(python)