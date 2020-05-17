project(Caribou.Geometry.Test)

set(SOURCE_FILES
        main.cpp)

enable_testing()

add_executable(${PROJECT_NAME} ${SOURCE_FILES})
target_link_libraries(${PROJECT_NAME} ${GTEST_BOTH_LIBRARIES})
target_link_libraries(${PROJECT_NAME} Caribou::Geometry)
