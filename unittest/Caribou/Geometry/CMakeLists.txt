project(Caribou.unittests.Geometry)

set(SOURCE_FILES
    main.cpp
)

if (NOT WIN32)
    find_package(Threads QUIET)
endif()

enable_testing()

add_executable(${PROJECT_NAME} ${SOURCE_FILES})
target_link_libraries(${PROJECT_NAME} gtest)
target_link_libraries(${PROJECT_NAME} Caribou::Geometry Caribou::Config)

if (NOT WIN32)
    target_link_libraries(${PROJECT_NAME} pthread)
elseif(MSVC)
    target_compile_options(${PROJECT_NAME} PRIVATE /bigobj)
endif()

list(APPEND target_rpath
    "$ORIGIN/../lib"
    "$ORIGIN/../../../lib"
    "@executable_path/../lib"
    "@executable_path/../../../lib"
)

set_target_properties(${PROJECT_NAME} PROPERTIES INSTALL_RPATH "${target_rpath}" )

install(
    TARGETS ${PROJECT_NAME}
    EXPORT Caribou
)