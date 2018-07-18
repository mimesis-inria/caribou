cmake_minimum_required(VERSION 3.1)
project(SofaModule)

set(LIB_NAME sofa)

set(HEADER_FILES
    macros.h
    helper.h
    core/__init__.h
    core/objetmodel/__init__.h
    core/objetmodel/Base.h
    simulation/__init__.h
    simulation/graph/__init__.h
    simulation/Simulation.h
)

set(SOURCE_FILES
    core/__init__.cpp
    core/objetmodel/__init__.cpp
    core/objetmodel/Base.cpp
    simulation/__init__.cpp
    simulation/graph/__init__.cpp
    sofa.cpp
    simulation/Base.cpp
)

find_package(SofaGeneral REQUIRED)
find_package(SofaBase REQUIRED)

# OPTIONS
find_package(PythonLibs 2 QUIET)
option(SOFA_PYTHON_WITH_PYTHON_2 "Compile the plugin with python 2 support" ${PYTHONLIBS_FOUND})

find_package(PythonLibs 3 QUIET)
option(SOFA_PYTHON_WITH_PYTHON_3 "Compile the plugin with python 3 support" ${PYTHONLIBS_FOUND})

if (NOT SOFA_PYTHON_WITH_PYTHON_2 AND NOT SOFA_PYTHON_WITH_PYTHON_3)
    message(SEND_ERROR No python package found. Either python 2 or python 3 must be provided for this module.)
else ()
    if (SOFA_PYTHON_WITH_PYTHON_2)
        message(STATUS " WITH PYTHON 2 SUPPORT (version ${PYTHONLIBS_VERSION_STRING})")
        find_package(PythonLibs 2 REQUIRED)
        include_directories(${PYTHON_INCLUDE_DIRS})
        set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/site-packages/2)
        add_definitions("-DPYTHON_2")
        add_library(${LIB_NAME} SHARED ${HEADER_FILES} ${SOURCE_FILES} ${SCENE_FILES})
        target_include_directories(${LIB_NAME}
            PRIVATE
            $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}>
            PUBLIC
            $<INSTALL_INTERFACE:include/caribou>
            )
        set_target_properties(
            ${LIB_NAME}
            PROPERTIES
            OUTPUT_NAME ${LIB_NAME}
            PREFIX ""
            SUFFIX ".so"
        )
        target_link_libraries(${LIB_NAME} SofaComponentGeneral SofaBaseTopology ${PYTHON_LIBRARIES})
    endif()

    if (SOFA_PYTHON_WITH_PYTHON_3)
        message(STATUS " WITH PYTHON 3 SUPPORT (version ${PYTHONLIBS_VERSION_STRING})")
        find_package(PythonLibs 3 REQUIRED)
        include_directories(${PYTHON_INCLUDE_DIRS})
        set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/site-packages/3)
        add_definitions("-DPYTHON_3")
        add_library(${LIB_NAME} SHARED ${HEADER_FILES} ${SOURCE_FILES} ${SCENE_FILES})
        set_target_properties(
            ${LIB_NAME}
            PROPERTIES
            OUTPUT_NAME ${LIB_NAME}
            PREFIX ""
            SUFFIX ".so"
        )
        target_link_libraries(${LIB_NAME} SofaComponentGeneral SofaBaseTopology ${PYTHON_LIBRARIES})
    endif()
endif ()

