cmake_minimum_required(VERSION 3.1)
project(SofaCaribou)

set(HEADER_FILES
    GraphComponents/Forcefield/HexahedronElasticForce.h
    GraphComponents/Forcefield/TetrahedronElasticForce.h
    GraphComponents/Forcefield/TractionForce.h
    GraphComponents/Ode/StaticODESolver.h
    GraphComponents/Topology/FictitiousGrid.h
    Traits.h)

set(TEMPLATE_FILES
    GraphComponents/Topology/FictitiousGrid.inl
        )

set(SOURCE_FILES
    GraphComponents/Forcefield/HexahedronElasticForce.cpp
    GraphComponents/Forcefield/TetrahedronElasticForce.cpp
    GraphComponents/Forcefield/TractionForce.cpp
    GraphComponents/Ode/StaticODESolver.cpp
    GraphComponents/Topology/FictitiousGrid.cpp
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
find_package(OpenMP)
find_package(LAPACK)
find_package(LAPACKE)
find_package(MKL)

add_library(${PROJECT_NAME} SHARED ${SOURCE_FILES})
target_include_directories(${PROJECT_NAME} PUBLIC "$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src/>")
target_include_directories(${PROJECT_NAME} PUBLIC "$<INSTALL_INTERFACE:include>")
target_include_directories(${PROJECT_NAME} SYSTEM PUBLIC Caribou.Topology Caribou.Geometry)
target_include_directories(${PROJECT_NAME} SYSTEM PUBLIC SofaCore)
target_link_libraries(${PROJECT_NAME} PUBLIC SofaCore SofaBaseTopology)
target_link_libraries(${PROJECT_NAME} PUBLIC Caribou.Topology Caribou.Geometry)


if (LAPACK_FOUND AND LAPACKE_FOUND AND MKL_FOUND)
    OPTION(Caribou_with_EIGEN_MKL "Use EIGEN with MKL/LAPACK/LAPACKE support" ON)
else()
    OPTION(Caribou_with_EIGEN_MKL "Use EIGEN with MKL/LAPACK/LAPACKE support" OFF)
endif()

if (OPENMP_FOUND)
    OPTION(Caribou_with_OpenMP "Use OpenMP support" ON)
else()
    OPTION(Caribou_with_OpenMP "Use OpenMP support" OFF)
endif()


if (Caribou_with_OpenMP)
    find_package(OpenMP REQUIRED)
    target_link_libraries(${PROJECT_NAME} PUBLIC ${OpenMP_CXX_LIBRARIES})
endif()

if (Caribou_with_EIGEN_MKL)
    find_package(LAPACK REQUIRED)
    find_package(LAPACKE REQUIRED)
    find_package(MKL REQUIRED)
    target_compile_definitions(${PROJECT_NAME} PUBLIC EIGEN_USE_LAPACKE)
    target_include_directories(${PROJECT_NAME} PUBLIC ${LAPACKE_INCLUDE_DIRS} ${LAPACKE_INCLUDE_DIRS_DEP} ${MKL_INCLUDE_DIRS})
    target_link_libraries(${PROJECT_NAME} PUBLIC ${LAPACK_LIBRARIES} ${LAPACKE_LIBRARIES} ${LAPACKE_LIBRARIES_DEP} ${MKL_LIBRARIES})
endif()


if (CARIBOU_WITH_PYTHON_2 OR CARIBOU_WITH_PYTHON_3)
    add_subdirectory(Python)
endif()