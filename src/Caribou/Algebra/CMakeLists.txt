cmake_minimum_required(VERSION 3.1)
project(Caribou.Algebra)

set(HEADER_FILES
        Vector.h
        Matrix.h
        SquareMatrix.h
        )

add_library(${PROJECT_NAME} INTERFACE)

target_include_directories(${PROJECT_NAME} INTERFACE "$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src/>")
target_include_directories(${PROJECT_NAME} INTERFACE "$<INSTALL_INTERFACE:include>")

if (CARIBOU_BUILD_TESTS)
    add_subdirectory(test)
endif()