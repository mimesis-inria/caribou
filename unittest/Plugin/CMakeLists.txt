project(SofaCaribouTest)

set(SOURCE_FILES
        main.cpp
        Algebra/test_base_vector_operations.cpp
        Algebra/test_eigen_matrix_wrapper.cpp
        Algebra/test_eigen_vector_wrapper.cpp
        Forcefield/test_tractionforce.cpp
        ODE/test_static.cpp
        Topology/test_fictitiousgrid.cpp
)

enable_testing()

find_package(SOFA COMPONENTS SofaFramework SofaSimulationGraph SofaBaseMechanics SofaBaseUtils SofaBaseLinearSolver QUIET REQUIRED)

add_executable(${PROJECT_NAME} ${SOURCE_FILES})
target_link_libraries(${PROJECT_NAME} PUBLIC gtest pthread)
target_link_libraries(${PROJECT_NAME} PUBLIC SofaCaribou)
target_link_libraries(${PROJECT_NAME} PUBLIC SofaHelper SofaSimulationGraph SofaBaseMechanics SofaBaseUtils SofaSimulationCommon SofaBaseLinearSolver)
target_link_libraries(${PROJECT_NAME} PUBLIC
                      $<$<AND:$<PLATFORM_ID:Linux>,$<VERSION_LESS:$<CXX_COMPILER_VERSION>,8.0>>:stdc++fs>
                      $<$<AND:$<PLATFORM_ID:Darwin>,$<VERSION_LESS:$<CXX_COMPILER_VERSION>,9.0>>:c++fs>)
target_compile_definitions(${PROJECT_NAME} PUBLIC
                           $<$<AND:$<PLATFORM_ID:Linux>,$<VERSION_LESS:$<CXX_COMPILER_VERSION>,8.0>>:LEGACY_CXX>
                           $<$<AND:$<PLATFORM_ID:Darwin>,$<VERSION_LESS:$<CXX_COMPILER_VERSION>,9.0>>:LEGACY_CXX>)

target_include_directories(${PROJECT_NAME} PUBLIC "$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src/>")

list(APPEND target_rpath
    "$ORIGIN/../lib"
    "@executable_path/../lib"
)

set_target_properties(${PROJECT_NAME} PROPERTIES INSTALL_RPATH "${target_rpath}" )

install(
        TARGETS ${PROJECT_NAME}
        EXPORT SofaCaribouTargets
        RUNTIME DESTINATION "bin"
)