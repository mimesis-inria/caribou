project(SofaCaribou VERSION 20.12.00)

set(HEADER_FILES
    config.h.in
    Algebra/BaseVectorOperations.h
    Algebra/EigenMatrix.h
    Algebra/EigenVector.h
    Forcefield/FictitiousGridElasticForce.h
    Forcefield/FictitiousGridHyperelasticForce.h
    Forcefield/HexahedronElasticForce.h
    Forcefield/HyperelasticForcefield.h
    Forcefield/TetrahedronElasticForce.h
    Forcefield/TractionForce.h
    Material/HyperelasticMaterial.h
    Material/NeoHookeanMaterial.h
    Material/SaintVenantKirchhoffMaterial.h
    Ode/BackwardEulerODESolver.h
    Ode/LegacyStaticODESolver.h
    Ode/NewtonRaphsonSolver.h
    Ode/StaticODESolver.h
    Solver/ConjugateGradientSolver.h
    Solver/EigenSolver.h
    Solver/LDLTSolver.h
    Solver/LinearSolver.h
    Solver/LLTSolver.h
    Solver/LUSolver.h
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
    Solver/EigenSolver.inl
    Topology/FictitiousGrid.inl
    Solver/LDLTSolver.inl
    Solver/LLTSolver.inl
    Solver/LUSolver.inl
    )

set(SOURCE_FILES
    Algebra/BaseVectorOperations.cpp
    Forcefield/FictitiousGridElasticForce.cpp
    Forcefield/FictitiousGridHyperelasticForce.cpp
    Forcefield/HexahedronElasticForce.cpp
    Forcefield/HyperelasticForcefield.cpp
    Forcefield/TetrahedronElasticForce.cpp
    Forcefield/TractionForce.cpp
    Material/HyperelasticMaterial.cpp
    Ode/BackwardEulerODESolver.cpp
    Ode/LegacyStaticODESolver.cpp
    Ode/NewtonRaphsonSolver.cpp
    Ode/StaticODESolver.cpp
    Solver/ConjugateGradientSolver.cpp
    Solver/LDLTSolver.cpp
    Solver/LLTSolver.cpp
    Solver/LUSolver.cpp
    Topology/FictitiousGrid.cpp
    Topology/IsoSurface.cpp
    Visitor/AssembleGlobalMatrix.cpp
    Visitor/ConstrainGlobalMatrix.cpp
    Visitor/MultiVecEqualVisitor.cpp
    init.cpp
)

# Caribou
if (NOT TARGET Caribou::Algebra)
    find_package(Caribou COMPONENTS Algebra Geometry Topology Mechanics REQUIRED QUIET)
endif()

# Sofa's packages
find_package(SOFA COMPONENTS SofaFramework SofaBaseLinearSolver SofaBaseTopology SofaEigen2Solver QUIET REQUIRED)

# Optional packages
find_package(Eigen3 QUIET REQUIRED)

add_library(${PROJECT_NAME} SHARED ${SOURCE_FILES} ${TEMPLATE_FILES} ${HEADER_FILES})
add_library(Caribou::Plugin ALIAS ${PROJECT_NAME})

target_link_libraries(${PROJECT_NAME} PUBLIC SofaCore SofaBaseTopology SofaBaseLinearSolver SofaEigen2Solver)
target_link_libraries(${PROJECT_NAME} PUBLIC Caribou::Algebra Caribou::Geometry Caribou::Topology Caribou::Mechanics)
target_link_libraries(${PROJECT_NAME} PUBLIC Eigen3::Eigen)
target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_17)
set_target_properties(${PROJECT_NAME} PROPERTIES POSITION_INDEPENDENT_CODE ON)

if (CARIBOU_OPTIMIZE_FOR_NATIVE)
#    todo(jnbrunet2000@gmail.com) Deactivated since it seems to create mis-alignments
#                                 issues for the SparseMatrix used in the CG. To reproduce:
#                                 https://gist.github.com/jnbrunet/7dfb354fbf2c6e0971674a1a83c29806
#    target_compile_options(${PROJECT_NAME} PRIVATE -march=native)
endif()

if (CARIBOU_WITH_MKL)
    set(MKL_STATIC ON)
    if (NOT CARIBOU_WITH_OPENMP)
        set(MKL_THREADING_VENDOR SEQUENTIAL)
    endif()
    find_package(MKL REQUIRED QUIET)

    target_compile_definitions(${PROJECT_NAME} PUBLIC CARIBOU_WITH_MKL)
    target_compile_definitions(${PROJECT_NAME} PRIVATE EIGEN_USE_MKL_ALL)
    target_include_directories(${PROJECT_NAME} PRIVATE ${MKL_INCLUDE_DIRS})
    target_link_libraries(${PROJECT_NAME} PRIVATE  ${MKL_LIBRARIES})
    target_link_options(${PROJECT_NAME} PRIVATE ${MKL_LINKER_FLAGS})
endif()


if (CARIBOU_WITH_OPENMP)
    find_package(OpenMP REQUIRED QUIET)
    target_link_libraries(${PROJECT_NAME} PUBLIC OpenMP::OpenMP_CXX)
    target_compile_definitions(${PROJECT_NAME} PUBLIC CARIBOU_WITH_OPENMP)
endif()

# Add the target to the component of the same name. This will enable to do:
# CMakeLists.txt
#    find_package(SofaCaribou)
#    add_library(my_cool_plugin SHARED)
#    target_link_libraries(my_cool_library PRIVATE SofaCaribou)
caribou_add_target_to_component (
    TARGET_NAME        ${PROJECT_NAME}
    COMPONENT_NAME     SofaCaribou
    PACKAGE_NAME       SofaCaribou
    PUBLIC_HEADERS     "${TEMPLATE_FILES};${HEADER_FILES}"
    HEADER_SRC_PREFIX  "${CMAKE_CURRENT_SOURCE_DIR}/.."
)

find_package(SofaPython3 QUIET)
if (CARIBOU_WITH_SOFA AND SofaPython3_FOUND)
    add_subdirectory(Python)
    set(CARIBOU_WITH_SP3 ON)
endif()

install (
    EXPORT SofaCaribouTargets
    DESTINATION "lib/cmake/SofaCaribou"
    COMPONENT headers
)

# This is required to be able to do "find_package(SofaCaribou)"
caribou_create_package(SofaCaribou ${CARIBOU_VERSION})
