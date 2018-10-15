cmake_minimum_required(VERSION 3.8)
project(Caribou.Components.Topology)

set(HEADER_FILES
    MultilevelSparseGridTopology.h
)

set(SOURCE_FILES
    MultilevelSparseGridTopology.cpp
)

find_package(SofaBase REQUIRED)

get_target_property(SOFA_CORE_INCLUDE_PATH SofaCore INTERFACE_INCLUDE_DIRECTORIES)

add_library(${PROJECT_NAME} SHARED ${SOURCE_FILES} ${HEADER_FILES})

target_include_directories(${PROJECT_NAME} SYSTEM PUBLIC ${SOFA_CORE_INCLUDE_PATH})
target_link_libraries(${PROJECT_NAME} SofaBaseTopology)

get_target_property(SOFA_CORE_INCLUDE_PATH SofaCore INTERFACE_INCLUDE_DIRECTORIES)
message(STATUS "SofaCore prefix: ${SOFA_CORE_INCLUDE_PATH}")