add_subdirectory(Config)
add_subdirectory(Geometry)
add_subdirectory(Topology)

find_package(SofaPython3 QUIET)
if (CARIBOU_WITH_SOFA AND SofaPython3_FOUND)
    add_subdirectory(Plugin)
endif()