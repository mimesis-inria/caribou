cmake_minimum_required(VERSION 3.1)
project(Algebra)

set(HEADER_FILES
        Vector.h
        )

set(SOURCE_FILES
        )

if (CARIBOU_BUILD_TESTS)
    add_subdirectory(test)
endif()