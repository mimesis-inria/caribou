project(SofaCaribou)

set(HEADER_FILES
    Algebra/EigenMatrixWrapper.h
    Forcefield/FictitiousGridElasticForce.h
    Forcefield/FictitiousGridHyperelasticForce.h
    Forcefield/HexahedronElasticForce.h
    Forcefield/HyperelasticForcefield.h
    Forcefield/TetrahedronElasticForce.h
    Forcefield/TractionForce.h
    Material/HyperelasticMaterial.h
    Material/NeoHookeanMaterial.h
    Material/SaintVenantKirchhoffMaterial.h
    Ode/StaticODESolver.h
    Solver/ConjugateGradientSolver.h
    Topology/CircleIsoSurface.h
    Topology/CylinderIsoSurface.h
    Topology/FictitiousGrid.h
    Topology/IsoSurface.h
    Topology/SphereIsoSurface.h
    Visitor/AssembleGlobalMatrix.h
    Visitor/ConstrainGlobalMatrix.h
    Visitor/MultiVecEqualVisitor.h
    )

set(TEMPLATE_FILES
    Forcefield/HyperelasticForcefield.inl
    Topology/FictitiousGrid.inl
    )

set(SOURCE_FILES
    Forcefield/FictitiousGridElasticForce.cpp
    Forcefield/FictitiousGridHyperelasticForce.cpp
    Forcefield/HexahedronElasticForce.cpp
    Forcefield/HyperelasticForcefield.cpp
    Forcefield/TetrahedronElasticForce.cpp
    Forcefield/TractionForce.cpp
    Material/HyperelasticMaterial.cpp
    Ode/StaticODESolver.cpp
    Solver/ConjugateGradientSolver.cpp
    Topology/FictitiousGrid.cpp
    Topology/IsoSurface.cpp
    Visitor/AssembleGlobalMatrix.cpp
    Visitor/ConstrainGlobalMatrix.cpp
    Visitor/MultiVecEqualVisitor.cpp
    init.cpp)

# Try to find Sofa's cmake prefix
# First, make sure it isn't already been manually set by the user's cmake variables
find_path (SOFA_PREFIX SofaFramework/SofaFrameworkConfig.cmake)
if (SOFA_PREFIX)
    find_package(SofaFramework QUIET)
    if (NOT SofaFramework_FOUND AND NOT "${SOFA_PREFIX}" IN_LIST CMAKE_PREFIX_PATH)
        list(APPEND CMAKE_PREFIX_PATH "${SOFA_PREFIX}")
    endif()
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

# Caribou
if (NOT TARGET Caribou::Algebra)
    find_package(Caribou COMPONENTS Algebra Geometry Topology Mechanics REQUIRED)
endif()

# Sofa's packages
find_package(SofaFramework REQUIRED)
find_package(SofaBase REQUIRED)
find_package(SofaCommon REQUIRED)

# Optional packages
find_package(OpenMP)
find_package(LAPACK)
find_package(LAPACKE)
find_package(MKL)
find_package(Eigen3 REQUIRED NO_MODULE)

add_library(${PROJECT_NAME} SHARED ${SOURCE_FILES} ${TEMPLATE_FILES} ${HEADER_FILES})
target_include_directories(${PROJECT_NAME} INTERFACE "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../../>")
target_include_directories(${PROJECT_NAME} INTERFACE "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/../../>")
target_include_directories(${PROJECT_NAME} PUBLIC "$<INSTALL_INTERFACE:include>")
target_link_libraries(${PROJECT_NAME} PUBLIC SofaCore SofaBaseTopology SofaBaseLinearSolver SofaEigen2Solver)
target_link_libraries(${PROJECT_NAME} PUBLIC Caribou::Algebra Caribou::Geometry Caribou::Topology Caribou::Mechanics)
set_property(TARGET ${PROJECT_NAME} PROPERTY CXX_STANDARD 17)
target_compile_features(${PROJECT_NAME} INTERFACE cxx_std_17)

if (CARIBOU_OPTIMIZE_FOR_NATIVE)
    target_compile_options(${PROJECT_NAME} PRIVATE)
endif()

if (LAPACK_FOUND AND LAPACKE_FOUND AND MKL_FOUND)
    OPTION(CARIBOU_WITH_EIGEN_MKL "Use EIGEN with MKL/LAPACK/LAPACKE support" ON)
else()
    OPTION(CARIBOU_WITH_EIGEN_MKL "Use EIGEN with MKL/LAPACK/LAPACKE support" OFF)
endif()

if (CARIBOU_WITH_EIGEN_MKL)
    find_package(LAPACK REQUIRED)
    find_package(LAPACKE REQUIRED)
    find_package(MKL REQUIRED)
    target_compile_definitions(${PROJECT_NAME} PUBLIC EIGEN_USE_LAPACKE)
    target_include_directories(${PROJECT_NAME} PUBLIC ${LAPACKE_INCLUDE_DIRS} ${LAPACKE_INCLUDE_DIRS_DEP} ${MKL_INCLUDE_DIRS})
    target_link_libraries(${PROJECT_NAME} PUBLIC ${LAPACK_LIBRARIES} ${LAPACKE_LIBRARIES} ${LAPACKE_LIBRARIES_DEP} ${MKL_LIBRARIES})
endif()


if (CARIBOU_WITH_OPENMP)
    find_package(OpenMP REQUIRED)
    target_link_libraries(${PROJECT_NAME} PUBLIC OpenMP::OpenMP_CXX)
    target_compile_definitions(${PROJECT_NAME} PUBLIC CARIBOU_WITH_OPENMP)
endif()

if (CARIBOU_BUILD_TESTS)
    add_subdirectory(test)
endif()

if (CARIBOU_WITH_PYTHON_2 OR CARIBOU_WITH_PYTHON_3)
    add_subdirectory(Python)
endif()

# CMake files
list(APPEND HEADER_FILES ${TEMPLATE_FILES})
foreach(header ${HEADER_FILES})
    file(RELATIVE_PATH path_from_package "${CMAKE_CURRENT_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}/${header}")
    get_filename_component(output_directory "${path_from_package}" DIRECTORY)
    if (header MATCHES "\\.in$")
        string(REGEX REPLACE "\\.in$" "" output_filepath ${path_from_package})
        configure_file("${header}" "${CMAKE_CURRENT_BINARY_DIR}/${output_filepath}")
        install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${output_filepath} DESTINATION "include/SofaCaribou/${output_directory}" COMPONENT headers)
    else()
        install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/${path_from_package} DESTINATION "include/SofaCaribou/${output_directory}" COMPONENT headers)
    endif()
endforeach()

install(TARGETS ${PROJECT_NAME}
    EXPORT ${PROJECT_NAME}Targets
    RUNTIME DESTINATION "bin" COMPONENT applications
    LIBRARY DESTINATION "lib" COMPONENT libraries
    ARCHIVE DESTINATION "lib" COMPONENT libraries
    PUBLIC_HEADER DESTINATION "include/${PROJECT_NAME}" COMPONENT headers
)

configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/SofaCaribouConfig.cmake.in"
    "${CMAKE_BINARY_DIR}/cmake/SofaCaribou/SofaCaribouConfig.cmake"
    INSTALL_DESTINATION
    lib/cmake/SofaCaribou
)
write_basic_package_version_file(
        ${CMAKE_BINARY_DIR}/cmake/SofaCaribou/SofaCaribouConfigVersion.cmake
        VERSION ${Caribou_VERSION}
        COMPATIBILITY AnyNewerVersion)

install(
    EXPORT SofaCaribouTargets 
    DESTINATION "lib/cmake/SofaCaribou" 
    COMPONENT headers)
install(FILES 
        "${CMAKE_BINARY_DIR}/cmake/SofaCaribou/SofaCaribouConfig.cmake"
        "${CMAKE_BINARY_DIR}/cmake/SofaCaribou/SofaCaribouConfigVersion.cmake"
        DESTINATION 
        "lib/cmake/SofaCaribou")
