cmake_minimum_required(VERSION 3.1)
project(Topology)

set(HEADER_FILES
        )

set(SOURCE_FILES
        )

if (CARIBOU_BUILD_TESTS)
    add_subdirectory(test)
endif()