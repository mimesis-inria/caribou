cmake_minimum_required(VERSION 3.1)
project(Geometry)

set(HEADER_FILES
        Entity.h
        Point.h
        )

set(SOURCE_FILES
        )

#add_library(${PROJECT_NAME} SHARED ${HEADER_FILES} ${SOURCE_FILES})

#target_include_directories(${PROJECT_NAME} PUBLIC "$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src/>")
#target_include_directories(${PROJECT_NAME} PUBLIC "$<INSTALL_INTERFACE:include>")

if (CARIBOU_BUILD_TESTS)
    add_subdirectory(test)
endif()

if (CARIBOU_WITH_PYTHON_2 OR CARIBOU_WITH_PYTHON_3)
    add_subdirectory(Python)
endif()