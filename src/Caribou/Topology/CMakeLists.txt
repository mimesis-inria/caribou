cmake_minimum_required(VERSION 3.1)
project(Topology)

set(HEADER_FILES
    Grid.h
)

set(TEMPLATE_FILES
    Grid.inl
)

set(SOURCE_FILES
    Grid.cpp
)

add_library(${PROJECT_NAME} SHARED ${HEADER_FILES} ${TEMPLATE_FILES} ${SOURCE_FILES})

target_include_directories(${PROJECT_NAME} PUBLIC "$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src/>")
target_include_directories(${PROJECT_NAME} PUBLIC "$<INSTALL_INTERFACE:include>")

if (CARIBOU_BUILD_TESTS)
    add_subdirectory(test)
endif()