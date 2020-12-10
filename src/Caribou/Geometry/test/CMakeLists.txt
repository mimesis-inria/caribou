project(Caribou.Geometry.Test)

set(SOURCE_FILES
        main.cpp)

enable_testing()

add_executable(${PROJECT_NAME} ${SOURCE_FILES})
target_link_libraries(${PROJECT_NAME} gtest pthread)
target_link_libraries(${PROJECT_NAME} Caribou::Geometry)

set_target_properties(${PROJECT_NAME} PROPERTIES INSTALL_RPATH "$ORIGIN/../lib" )

install(
    TARGETS ${PROJECT_NAME}
    EXPORT Caribou
)