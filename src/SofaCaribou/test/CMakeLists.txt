project(SofaCaribouTest)

set(SOURCE_FILES
        main.cpp
        Algebra/test_eigen_matrix_wrapper.cpp
        Forcefield/test_tractionforce.cpp
)

enable_testing()

find_package(SofaFramework QUIET REQUIRED)
find_package(SofaSimulation QUIET REQUIRED)
find_package(SofaBase QUIET REQUIRED)

add_executable(${PROJECT_NAME} ${SOURCE_FILES})
target_link_libraries(${PROJECT_NAME} ${GTEST_BOTH_LIBRARIES})
target_link_libraries(${PROJECT_NAME} SofaCaribou)
target_link_libraries(${PROJECT_NAME} SofaHelper SofaSimulationGraph SofaBaseMechanics)

target_include_directories(${PROJECT_NAME} PUBLIC "$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src/>")
target_include_directories(${PROJECT_NAME} PUBLIC "$<BUILD_INTERFACE:${GTEST_INCLUDE_DIR}>")
