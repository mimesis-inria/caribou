add_subdirectory(Algebra)
add_subdirectory(Config)
add_subdirectory(Geometry)
add_subdirectory(Mechanics)
add_subdirectory(Topology)

if (CARIBOU_WITH_PYTHON_3)
    add_subdirectory(Python)
endif()

caribou_create_package(Caribou ${CARIBOU_VERSION})