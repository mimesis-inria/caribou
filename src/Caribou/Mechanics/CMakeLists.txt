project(Mechanics)

set(HEADER_FILES
        Elasticity/Strain.h)

find_package(Eigen3 QUIET REQUIRED NO_MODULE)

add_library(${PROJECT_NAME} INTERFACE)
add_library(Caribou::${PROJECT_NAME} ALIAS ${PROJECT_NAME})

target_include_directories(${PROJECT_NAME} INTERFACE "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../../>")
target_include_directories(${PROJECT_NAME} INTERFACE "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/../../>")
target_include_directories(${PROJECT_NAME} INTERFACE "$<INSTALL_INTERFACE:include>")
target_link_libraries(${PROJECT_NAME} INTERFACE Eigen3::Eigen)
target_compile_features(${PROJECT_NAME} INTERFACE cxx_std_14)

caribou_install_target(Caribou ${PROJECT_NAME} ${HEADER_FILES})

if (CARIBOU_BUILD_TESTS)
    add_subdirectory(test)
endif()
