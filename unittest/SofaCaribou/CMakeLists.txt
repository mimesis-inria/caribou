project(Caribou.unittests.SofaCaribou)

set(SOURCE_FILES
        main.cpp
        Algebra/test_base_vector_operations.cpp
        Algebra/test_eigen_matrix_wrapper.cpp
        Algebra/test_eigen_vector_wrapper.cpp
        Forcefield/test_hyperelasticforcefield.cpp
        Forcefield/test_tractionforce.cpp
        Mass/test_cariboumass.cpp
        ODE/test_backward_euler.cpp
        ODE/test_static.cpp
        Topology/test_fictitiousgrid.cpp
)

enable_testing()

find_package(SOFA COMPONENTS Sofa.Framework Sofa.Simulation.Graph Sofa.Component.StateContainer Sofa.Component.Mass Sofa.Component.SceneUtility Sofa.Component.LinearSolver.Iterative Sofa.Component.Mass Sofa.Component.SolidMechanics.FEM.Elastic QUIET REQUIRED)
find_package(Caribou COMPONENTS Topology REQUIRED)
find_package(Threads REQUIRED)

if (CARIBOU_WITH_VTK)
    list(APPEND SOURCE_FILES Topology/test_cariboutopology.cpp)
endif()

add_executable(${PROJECT_NAME} ${SOURCE_FILES})

if (NOT WIN32)
    target_link_libraries(${PROJECT_NAME} PUBLIC pthread)
endif()

target_link_libraries(${PROJECT_NAME} PUBLIC gtest)
target_link_libraries(${PROJECT_NAME} PUBLIC SofaCaribou)
target_link_libraries(${PROJECT_NAME} PUBLIC Sofa.Helper Sofa.Component.StateContainer Sofa.Component.Mass Sofa.Component.SceneUtility Sofa.Component.LinearSolver.Iterative Sofa.Component.SolidMechanics.FEM.Elastic)
target_link_libraries(${PROJECT_NAME} PUBLIC
                      $<$<AND:$<PLATFORM_ID:Linux>,$<VERSION_LESS:$<CXX_COMPILER_VERSION>,8.0>>:stdc++fs>
                      $<$<AND:$<PLATFORM_ID:Darwin>,$<VERSION_LESS:$<CXX_COMPILER_VERSION>,9.0>>:c++fs>)
if (${Sofa.Framework_VERSION} VERSION_LESS "22.06.99")
    target_link_libraries(${PROJECT_NAME} PUBLIC SofaSimulationGraph)
else()
    target_link_libraries(${PROJECT_NAME} PUBLIC Sofa.Simulation.Graph)
endif()
target_compile_definitions(${PROJECT_NAME} PUBLIC
                           $<$<AND:$<PLATFORM_ID:Linux>,$<VERSION_LESS:$<CXX_COMPILER_VERSION>,8.0>>:LEGACY_CXX>
                           $<$<AND:$<PLATFORM_ID:Darwin>,$<VERSION_LESS:$<CXX_COMPILER_VERSION>,9.0>>:LEGACY_CXX>)

target_include_directories(${PROJECT_NAME} PUBLIC "$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src/>")

if (SOFA_VERSION VERSION_GREATER_EQUAL "20.12.99")
    target_link_libraries(${PROJECT_NAME} PUBLIC Sofa.Testing)
endif()

list(APPEND target_rpath
    "$ORIGIN/../lib"
    "$ORIGIN/../../../lib"
    "@executable_path/../lib"
    "@executable_path/../../../lib"
    # SofaMiscForceField
    "$ORIGIN/../../../plugins/SofaMiscForceField/lib"
    "@executable_path/../../../plugins/SofaMiscForceField/lib"
    # SofaMiscTopology
    "$ORIGIN/../../../plugins/SofaMiscTopology/lib"
    "@executable_path/../../../plugins/SofaMiscTopology/lib"
    # SofaDeformable
    "$ORIGIN/../../../plugins/SofaDeformable/lib"
    "@executable_path/../../../plugins/SofaDeformable/lib"
)

set_target_properties(${PROJECT_NAME} PROPERTIES INSTALL_RPATH "${target_rpath}" )

install(
        TARGETS ${PROJECT_NAME}
        EXPORT SofaCaribouTargets
        RUNTIME DESTINATION "bin"
)