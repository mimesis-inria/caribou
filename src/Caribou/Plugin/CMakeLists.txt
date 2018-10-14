cmake_minimum_required(VERSION 3.1)
project(Caribou.Plugin)

set(SOURCE_FILES
    init.cpp
)

find_package(SofaFramework REQUIRED)

add_library(${PROJECT_NAME} SHARED ${SOURCE_FILES})
target_include_directories(${PROJECT_NAME} PUBLIC "$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src/>")
target_include_directories(${PROJECT_NAME} PUBLIC "$<INSTALL_INTERFACE:include>")

target_link_libraries(${PROJECT_NAME} Caribou.Components.Topology)
target_link_libraries(${PROJECT_NAME} SofaCore)