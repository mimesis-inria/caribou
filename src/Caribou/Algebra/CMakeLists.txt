project(Algebra)

set(HEADER_FILES
        Tensor.h)

find_package(Eigen3 REQUIRED NO_MODULE)

add_library(${PROJECT_NAME} INTERFACE)
add_library(Caribou::Algebra ALIAS ${PROJECT_NAME})

target_link_libraries(${PROJECT_NAME} INTERFACE Eigen3::Eigen)
target_include_directories(${PROJECT_NAME} INTERFACE "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../../>")
target_include_directories(${PROJECT_NAME} INTERFACE "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/../../>")
target_include_directories(${PROJECT_NAME} INTERFACE "$<INSTALL_INTERFACE:include>")
target_compile_features(${PROJECT_NAME} INTERFACE cxx_std_14)

caribou_install_target(Caribou ${PROJECT_NAME} ${HEADER_FILES})

