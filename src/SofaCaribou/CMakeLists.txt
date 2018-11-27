cmake_minimum_required(VERSION 3.1)
project(SofaCaribou)

set(HEADER_FILES
    GraphComponents/Topology/MultilevelSparseGridTopology.h)

set(TEMPLATE_FILES
    GraphComponents/Topology/MultilevelSparseGridTopology.inl)

set(SOURCE_FILES
    GraphComponents/Topology/MultilevelSparseGridTopology.cpp
    init.cpp)

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


find_package(SofaBase REQUIRED)
find_package(SofaFramework REQUIRED)

add_library(${PROJECT_NAME} SHARED ${SOURCE_FILES})
target_include_directories(${PROJECT_NAME} PUBLIC "$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src/>")
target_include_directories(${PROJECT_NAME} PUBLIC "$<INSTALL_INTERFACE:include>")
target_include_directories(${PROJECT_NAME} SYSTEM PUBLIC Caribou.Topology)
target_include_directories(${PROJECT_NAME} SYSTEM PUBLIC SofaCore)
target_link_libraries(${PROJECT_NAME} SofaCore SofaBaseTopology)
target_link_libraries(${PROJECT_NAME} Caribou.Topology)