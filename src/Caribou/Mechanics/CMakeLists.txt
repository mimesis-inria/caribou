cmake_minimum_required(VERSION 3.1)
project(Caribou.Mechanics)

set(HEADER_FILES
        Elasticity/Strain.h)

find_package(Eigen3 REQUIRED NO_MODULE)

add_library(${PROJECT_NAME} INTERFACE)

target_link_libraries(${PROJECT_NAME} INTERFACE Eigen3::Eigen)

target_include_directories(${PROJECT_NAME} INTERFACE "$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src/>")
target_include_directories(${PROJECT_NAME} INTERFACE "$<INSTALL_INTERFACE:include>")

if (CARIBOU_BUILD_TESTS)
    add_subdirectory(test)
endif()