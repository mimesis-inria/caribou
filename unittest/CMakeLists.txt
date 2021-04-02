add_subdirectory(Geometry)
add_subdirectory(Mechanics)
add_subdirectory(Topology)


if (CARIBOU_WITH_SOFA)
    find_package(Sofa.Testing QUIET)
    if (TARGET Sofa.Testing OR SOFAHELPER_HAVE_GTEST)
        add_subdirectory(Plugin)
    endif()
endif()