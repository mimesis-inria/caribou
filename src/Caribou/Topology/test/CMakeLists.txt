project(Caribou.Topology.Test)

set(SOURCE_FILES
        main.cpp)

enable_testing()
find_package(GTest REQUIRED)

add_executable(${PROJECT_NAME} ${SOURCE_FILES})
target_link_libraries(${PROJECT_NAME} ${GTEST_BOTH_LIBRARIES} Topology)
target_include_directories(${PROJECT_NAME} PUBLIC "$<INSTALL_INTERFACE:include>")
target_include_directories(${PROJECT_NAME} PUBLIC "$<BUILD_INTERFACE:${GTEST_INCLUDE_DIR}>")
