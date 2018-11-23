cmake_minimum_required(VERSION 3.1)

project(Caribou.Geometry.Test)

set(SOURCE_FILES
        main.cpp)

enable_testing()
find_package(GTest REQUIRED)

add_executable(${PROJECT_NAME} ${SOURCE_FILES})
target_link_libraries(${PROJECT_NAME} ${GTEST_BOTH_LIBRARIES})
target_link_libraries(${PROJECT_NAME} Caribou.Algebra)

target_include_directories(${PROJECT_NAME} PUBLIC "$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src/>")
target_include_directories(${PROJECT_NAME} PUBLIC "$<INSTALL_INTERFACE:include>")

target_include_directories(${PROJECT_NAME} PUBLIC "$<BUILD_INTERFACE:${GTEST_INCLUDE_DIR}>")