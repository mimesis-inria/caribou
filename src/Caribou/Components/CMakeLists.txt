cmake_minimum_required(VERSION 3.1)
project(Caribou.Components)

add_subdirectory(Topology)

#add_library(${PROJECT_NAME} SHARED $<TARGET_OBJECTS:Caribou.Components.Topology>)
#target_link_libraries(${PROJECT_NAME} PUBLIC Caribou.Components.Topology)