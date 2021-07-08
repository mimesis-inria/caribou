add_subdirectory(Caribou/Geometry)
add_subdirectory(Caribou/Mechanics)
add_subdirectory(Caribou/Topology)


if (CARIBOU_WITH_SOFA)
    find_package(Sofa.Testing QUIET)
    if (TARGET Sofa.Testing OR SOFAHELPER_HAVE_GTEST)
        add_subdirectory(SofaCaribou)
    endif()
endif()
