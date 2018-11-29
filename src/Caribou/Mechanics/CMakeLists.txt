cmake_minimum_required(VERSION 3.1)
project(Caribou.Mechanics)

set(HEADER_FILES
        Elasticity/Strain.h)

set(SOURCE_FILES
        )

if (CARIBOU_BUILD_TESTS)
    add_subdirectory(test)
endif()