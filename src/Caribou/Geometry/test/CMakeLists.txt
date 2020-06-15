project(Caribou.Geometry.Test)

set(SOURCE_FILES
        main.cpp)

enable_testing()

add_executable(${PROJECT_NAME} ${SOURCE_FILES})
target_link_libraries(${PROJECT_NAME} ${GTEST_BOTH_LIBRARIES})
target_include_directories(${PROJECT_NAME} SYSTEM PRIVATE ${GTEST_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} Caribou::Geometry)
