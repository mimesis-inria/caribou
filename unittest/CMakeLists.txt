add_subdirectory(Geometry)
add_subdirectory(Mechanics)
add_subdirectory(Topology)

if (CARIBOU_WITH_SOFA AND SOFAHELPER_HAVE_GTEST)
    add_subdirectory(Plugin)
endif()