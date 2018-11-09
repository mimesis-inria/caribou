cmake_minimum_required(VERSION 3.1)
project(Caribou.Topology)

set(HEADER_FILES
        Engine/Grid/Grid.h
)

set(TEMPLATE_FILES
        Engine/Grid/Grid.inl
)

set(SOURCE_FILES
        Engine/Grid/Grid2D.cpp
        Engine/Grid/Grid3D.cpp
)

add_library(${PROJECT_NAME} SHARED ${HEADER_FILES} ${TEMPLATE_FILES} ${SOURCE_FILES})

target_include_directories(${PROJECT_NAME} PUBLIC "$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src/>")
target_include_directories(${PROJECT_NAME} PUBLIC "$<INSTALL_INTERFACE:include>")

if (CARIBOU_BUILD_TESTS)
    add_subdirectory(test)
endif()