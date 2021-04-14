project(Caribou.unittests.Geometry)

set(SOURCE_FILES
    main.cpp
)

enable_testing()

add_executable(${PROJECT_NAME} ${SOURCE_FILES})
target_link_libraries(${PROJECT_NAME} gtest pthread)
target_link_libraries(${PROJECT_NAME} Caribou::Geometry Caribou::Config)

list(APPEND target_rpath
    "$ORIGIN/../lib"
    "@executable_path/../lib"
)

set_target_properties(${PROJECT_NAME} PROPERTIES INSTALL_RPATH "${target_rpath}" )

install(
    TARGETS ${PROJECT_NAME}
    EXPORT Caribou
)