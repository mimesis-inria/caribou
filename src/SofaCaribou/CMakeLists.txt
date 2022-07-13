project(SofaCaribou)

set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/FEniCS/Material)

set(HEADER_FILES
    config.h.in
    Algebra/BaseVectorOperations.h
    Algebra/EigenMatrix.h
    Algebra/EigenVector.h
    Forcefield/CaribouForcefield.h
    Forcefield/CaribouForcefield[Hexahedron].h
    Forcefield/CaribouForcefield[Hexahedron20].h
    Forcefield/CaribouForcefield[Hexahedron_FEniCS].h
    Forcefield/CaribouForcefield[Hexahedron_FEniCS20].h
    Forcefield/CaribouForcefield[Quad].h
    Forcefield/CaribouForcefield[Quad8].h
    Forcefield/CaribouForcefield[Tetrahedron].h
    Forcefield/CaribouForcefield[Tetrahedron10].h
    Forcefield/CaribouForcefield[Triangle].h
    Forcefield/CaribouForcefield[Triangle6].h
    Forcefield/HexahedronElasticForce.h
    Forcefield/HyperelasticForcefield.h
    Forcefield/HyperelasticForcefield[Hexahedron].h
    Forcefield/HyperelasticForcefield[Hexahedron20].h
    Forcefield/HyperelasticForcefield[Tetrahedron].h
    Forcefield/HyperelasticForcefield[Tetrahedron10].h
    Forcefield/HyperelasticForcefield[Triangle].h
    FEniCS/Forcefield/HyperelasticForcefield_FEniCS.h
    FEniCS/Forcefield/HyperelasticForcefield_FEniCS[Hexahedron_FEniCS].h
    FEniCS/Forcefield/HyperelasticForcefield_FEniCS[Hexahedron_FEniCS20].h
    FEniCS/Forcefield/HyperelasticForcefield_FEniCS[Tetrahedron].h
    FEniCS/Forcefield/HyperelasticForcefield_FEniCS[Tetrahedron10].h
    Forcefield/TetrahedronElasticForce.h
    Forcefield/TractionForcefield.h
    Forcefield/TractionForcefield[Quad].h
    Forcefield/TractionForcefield[Quad8].h
    Forcefield/TractionForcefield[Triangle].h
    Forcefield/TractionForcefield[Triangle6].h
    Mapping/CaribouBarycentricMapping.h
    Mapping/CaribouBarycentricMapping[Hexahedron].h
    Mapping/CaribouBarycentricMapping[Hexahedron20].h
    Mapping/CaribouBarycentricMapping[Tetrahedron].h
    Mapping/CaribouBarycentricMapping[Tetrahedron10].h
    Mass/CaribouMass.h
    Mass/CaribouMass[Tetrahedron].h
    Mass/CaribouMass[Tetrahedron10].h
    Mass/CaribouMass[Hexahedron].h
    Mass/CaribouMass[Hexahedron_FEniCS].h
    Mass/CaribouMass[Hexahedron20].h
    Mass/CaribouMass[Hexahedron_FEniCS20].h
    Material/HyperelasticMaterial.h
    Material/NeoHookeanMaterial.h
    Material/SaintVenantKirchhoffMaterial.h
    # FEniCS Material
    FEniCS/Material/FEniCS_Material.h
    FEniCS/Material/FEniCS_Material[Hexahedron_FEniCS].h
    FEniCS/Material/FEniCS_Material[Hexahedron_FEniCS20].h
    FEniCS/Material/FEniCS_Material[Tetrahedron].h
    FEniCS/Material/FEniCS_Material[Tetrahedron10].h
    # StVk FEniCS Generated code
    FEniCS/Material/FEniCS_Generated_code/SaintVenantKirchhoff_Tetra.h
    FEniCS/Material/FEniCS_Generated_code/SaintVenantKirchhoff_Tetra_Order2.h
    FEniCS/Material/FEniCS_Generated_code/SaintVenantKirchhoff_Hexa.h
    FEniCS/Material/FEniCS_Generated_code/SaintVenantKirchhoff_Hexa_Order2.h
    # NeoHooke FEniCS Generated code
    FEniCS/Material/FEniCS_Generated_code/NeoHooke_Tetra.h
    FEniCS/Material/FEniCS_Generated_code/NeoHooke_Tetra_Order2.h
    FEniCS/Material/FEniCS_Generated_code/NeoHooke_Hexa.h
    FEniCS/Material/FEniCS_Generated_code/NeoHooke_Hexa_Order2.h
    # MooneyRivlin FEniCS Generated code
    FEniCS/Material/FEniCS_Generated_code/MooneyRivlin_Tetra.h
    FEniCS/Material/FEniCS_Generated_code/MooneyRivlin_Tetra_Order2.h
    FEniCS/Material/FEniCS_Generated_code/MooneyRivlin_Hexa.h
    FEniCS/Material/FEniCS_Generated_code/MooneyRivlin_Hexa_Order2.h
    # Ogden FEniCS Generated code
    FEniCS/Material/FEniCS_Generated_code/Ogden_Tetra.h
    FEniCS/Material/FEniCS_Generated_code/Ogden_Tetra_Order2.h
    FEniCS/Material/FEniCS_Generated_code/Ogden_Hexa.h
    FEniCS/Material/FEniCS_Generated_code/Ogden_Hexa_Order2.h
    ##################################
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
    Topology/CaribouTopology.h
    Topology/CaribouTopology[Hexahedron].h
    Topology/CaribouTopology[Hexahedron20].h
    Topology/CaribouTopology[Hexahedron_FEniCS].h
    Topology/CaribouTopology[Hexahedron_FEniCS20].h
    Topology/CaribouTopology[Quad].h
    Topology/CaribouTopology[Quad8].h
    Topology/CaribouTopology[Tetrahedron].h
    Topology/CaribouTopology[Tetrahedron10].h
    Topology/CaribouTopology[Triangle].h
    Topology/CaribouTopology[Triangle6].h
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
    Forcefield/CaribouForcefield.inl
    Forcefield/HyperelasticForcefield.inl
    FEniCS/Forcefield/HyperelasticForcefield_FEniCS.inl
    Forcefield/TractionForcefield.inl
    Mapping/CaribouBarycentricMapping.inl
    Mass/CaribouMass.inl
    Solver/ConjugateGradientSolver.inl
    Solver/EigenSolver.inl
    Topology/CaribouTopology.inl
    Topology/FictitiousGrid.inl
    Solver/LDLTSolver.inl
    Solver/LLTSolver.inl
    Solver/LUSolver.inl
    )

set(SOURCE_FILES
    Algebra/BaseVectorOperations.cpp
    Forcefield/CaribouForcefield[Hexahedron].cpp
    Forcefield/CaribouForcefield[Hexahedron20].cpp
    Forcefield/CaribouForcefield[Hexahedron_FEniCS].cpp
    Forcefield/CaribouForcefield[Hexahedron_FEniCS20].cpp
    Forcefield/CaribouForcefield[Quad].cpp
    Forcefield/CaribouForcefield[Quad8].cpp
    Forcefield/CaribouForcefield[Tetrahedron].cpp
    Forcefield/CaribouForcefield[Tetrahedron10].cpp
    Forcefield/CaribouForcefield[Triangle].cpp
    Forcefield/CaribouForcefield[Triangle6].cpp
    Forcefield/HexahedronElasticForce.cpp
    Forcefield/HyperelasticForcefield[Hexahedron].cpp
    Forcefield/HyperelasticForcefield[Hexahedron20].cpp
    Forcefield/HyperelasticForcefield[Tetrahedron].cpp
    Forcefield/HyperelasticForcefield[Tetrahedron10].cpp
    Forcefield/HyperelasticForcefield[Triangle].cpp
    FEniCS/Forcefield/HyperelasticForcefield_FEniCS[Hexahedron_FEniCS].cpp
    FEniCS/Forcefield/HyperelasticForcefield_FEniCS[Hexahedron_FEniCS20].cpp
    FEniCS/Forcefield/HyperelasticForcefield_FEniCS[Tetrahedron].cpp
    FEniCS/Forcefield/HyperelasticForcefield_FEniCS[Tetrahedron10].cpp
    Forcefield/TetrahedronElasticForce.cpp
    Forcefield/TractionForcefield[Quad].cpp
    Forcefield/TractionForcefield[Quad8].cpp
    Forcefield/TractionForcefield[Triangle].cpp
    Forcefield/TractionForcefield[Triangle6].cpp
    Mapping/CaribouBarycentricMapping[Hexahedron].cpp
    Mapping/CaribouBarycentricMapping[Hexahedron20].cpp
    Mapping/CaribouBarycentricMapping[Tetrahedron].cpp
    Mapping/CaribouBarycentricMapping[Tetrahedron10].cpp
    Mass/CaribouMass[Tetrahedron].cpp
    Mass/CaribouMass[Tetrahedron10].cpp
    Mass/CaribouMass[Hexahedron].cpp
    Mass/CaribouMass[Hexahedron_FEniCS].cpp
    Mass/CaribouMass[Hexahedron20].cpp
    Mass/CaribouMass[Hexahedron_FEniCS20].cpp
    Material/HyperelasticMaterial.cpp
    # FEniCS Material
    FEniCS/Material/FEniCS_Material[Hexahedron_FEniCS].cpp
    FEniCS/Material/FEniCS_Material[Hexahedron_FEniCS20].cpp
    FEniCS/Material/FEniCS_Material[Tetrahedron].cpp
    FEniCS/Material/FEniCS_Material[Tetrahedron10].cpp
    # StVk FEniCS Generated code
    FEniCS/Material/FEniCS_Generated_code/SaintVenantKirchhoff_Tetra.c
    FEniCS/Material/FEniCS_Generated_code/SaintVenantKirchhoff_Tetra_Order2.c
    FEniCS/Material/FEniCS_Generated_code/SaintVenantKirchhoff_Hexa.c
    FEniCS/Material/FEniCS_Generated_code/SaintVenantKirchhoff_Hexa_Order2.c
    # NeoHooke FEniCS Generated code
    FEniCS/Material/FEniCS_Generated_code/NeoHooke_Tetra.c
    FEniCS/Material/FEniCS_Generated_code/NeoHooke_Tetra_Order2.c
    FEniCS/Material/FEniCS_Generated_code/NeoHooke_Hexa.c
    FEniCS/Material/FEniCS_Generated_code/NeoHooke_Hexa_Order2.c
    # MooneyRivlin FEniCS Generated code
    FEniCS/Material/FEniCS_Generated_code/MooneyRivlin_Tetra.c
    FEniCS/Material/FEniCS_Generated_code/MooneyRivlin_Tetra_Order2.c
    FEniCS/Material/FEniCS_Generated_code/MooneyRivlin_Hexa.c
    FEniCS/Material/FEniCS_Generated_code/MooneyRivlin_Hexa_Order2.c
    # Ogden FEniCS Generated code
    FEniCS/Material/FEniCS_Generated_code/Ogden_Tetra.c
    FEniCS/Material/FEniCS_Generated_code/Ogden_Tetra_Order2.c
    FEniCS/Material/FEniCS_Generated_code/Ogden_Hexa.c
    FEniCS/Material/FEniCS_Generated_code/Ogden_Hexa_Order2.c
    ##################################
    Ode/BackwardEulerODESolver.cpp
    Ode/LegacyStaticODESolver.cpp
    Ode/NewtonRaphsonSolver.cpp
    Ode/StaticODESolver.cpp
    Solver/ConjugateGradientSolver.cpp
    Solver/EigenSolver.cpp
    Solver/LDLTSolver.cpp
    Solver/LLTSolver.cpp
    Solver/LUSolver.cpp
    Topology/CaribouTopology[Hexahedron].cpp
    Topology/CaribouTopology[Hexahedron20].cpp
    Topology/CaribouTopology[Hexahedron_FEniCS].cpp
    Topology/CaribouTopology[Hexahedron_FEniCS20].cpp
    Topology/CaribouTopology[Quad].cpp
    Topology/CaribouTopology[Quad8].cpp
    Topology/CaribouTopology[Tetrahedron].cpp
    Topology/CaribouTopology[Tetrahedron10].cpp
    Topology/CaribouTopology[Triangle].cpp
    Topology/CaribouTopology[Triangle6].cpp
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

# Bugged dependencies
if (macos)
    find_package(Boost COMPONENTS system QUIET) # Bug introduced in SOFA v21.06 where this dependencies isn't propagated
endif()

add_library(${PROJECT_NAME} SHARED ${SOURCE_FILES} ${TEMPLATE_FILES} ${HEADER_FILES})
add_library(Caribou::Plugin ALIAS ${PROJECT_NAME})

target_link_libraries(${PROJECT_NAME} PUBLIC SofaCore SofaBaseTopology SofaBaseLinearSolver SofaEigen2Solver)
target_link_libraries(${PROJECT_NAME} PUBLIC Caribou::Algebra Caribou::Geometry Caribou::Topology Caribou::Mechanics)
target_link_libraries(${PROJECT_NAME} PUBLIC Eigen3::Eigen)
target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_17)
set_target_properties(${PROJECT_NAME} PROPERTIES POSITION_INDEPENDENT_CODE ON)

list(APPEND target_rpath
    "$ORIGIN/../lib"
    "@executable_path/../lib"
)

set_target_properties(${PROJECT_NAME} PROPERTIES INSTALL_RPATH "${target_rpath}" )

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
    target_link_libraries(${PROJECT_NAME} PRIVATE  ${MKL_LIBRARIES} -ldl)
    target_link_options(${PROJECT_NAME} PRIVATE ${MKL_LINKER_FLAGS})
endif()


if (CARIBOU_WITH_OPENMP)
    find_package(OpenMP REQUIRED QUIET)
    target_link_libraries(${PROJECT_NAME} PUBLIC OpenMP::OpenMP_CXX)
    target_compile_definitions(${PROJECT_NAME} PUBLIC CARIBOU_WITH_OPENMP)
endif()

set_target_properties(${PROJECT_NAME} PROPERTIES DEBUG_POSTFIX "_d")

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
if (CARIBOU_WITH_SOFA AND SofaPython3_FOUND AND CARIBOU_WITH_PYTHON_3)
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
# Make the target relocatable
set_target_properties(${PROJECT_NAME} PROPERTIES RELOCATABLE_INSTALL_DIR "plugins/SofaCaribou")
set_target_properties(${PROJECT_NAME} PROPERTIES EXPORT_PROPERTIES "RELOCATABLE_INSTALL_DIR")

if (SOFA_VERSION VERSION_GREATER_EQUAL "20.12.00")
    # Set RPATH towards relocatable dependencies
    sofa_auto_set_target_rpath(TARGETS ${PROJECT_NAME} RELOCATABLE "plugins")
endif()
