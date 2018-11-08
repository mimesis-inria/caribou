cmake_minimum_required(VERSION 3.1)
project(Caribou.Geometry.Python)

set(DIRECTORY Caribou/Geometry)

set(PYTHON_FILES
    __init__.py
)

set(PYTHON_TEST_FILES
    tests/test_point.py
)

set(HEADER_FILES
    Point.h
)

set(SOURCE_FILES
    Geometry.cpp
    Point.cpp
)

if (CARIBOU_WITH_PYTHON_2)
    find_package(PythonLibs 2.7 REQUIRED)
    message(STATUS "${PROJECT_NAME} with python 2.7 support (version ${PYTHONLIBS_VERSION_STRING})")

    set(PYBIND11_PYTHON_VERSION 2.7 CACHE STRING "")
    find_package(pybind11 REQUIRED)

    add_library(${PROJECT_NAME} SHARED ${SOURCE_FILES} ${PYTHON_FILES})

    target_include_directories(${PROJECT_NAME}
            PUBLIC ${PYTHON_INCLUDE_DIRS}
            PUBLIC "$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src/>"
            PUBLIC $<INSTALL_INTERFACE:include>
            )


    target_compile_definitions(${PROJECT_NAME}
            PUBLIC PYTHON_27
    )

    set_target_properties(
        ${PROJECT_NAME}
        PROPERTIES
            OUTPUT_NAME ${PROJECT_NAME}
            PREFIX "${PYTHON_MODULE_PREFIX}"
            SUFFIX "${PYTHON_MODULE_EXTENSION}"
            LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/site-packages/2.7/${DIRECTORY}"
    )

    target_link_libraries(${PROJECT_NAME} pybind11::module)

    set(PACKAGE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/site-packages/2.7")

    file(COPY ${PYTHON_FILES} DESTINATION "${PACKAGE_INSTALL_PREFIX}/${DIRECTORY}")
        foreach(t ${PYTHON_TEST_FILES})
            configure_file(${t} "${CMAKE_BINARY_DIR}/python_tests/2.7/${DIRECTORY}/${t}")
        endforeach()
    install(FILES ${PYTHON_FILES} DESTINATION site-packages/2.7/${DIRECTORY})
endif()