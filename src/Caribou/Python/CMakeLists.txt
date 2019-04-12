cmake_minimum_required(VERSION 3.1)
project(Caribou.Python)

set(DIRECTORY Caribou)

set(PYTHON_FILES
        __init__.py
        )

set(PYTHON_TEST_FILES
        )

set(SOURCE_FILES
        )

if (CARIBOU_WITH_PYTHON_2)
    caribou_add_python_module(CaribouPython
                              PYTHON_VERSION 2.7
                              DESTINATION Caribou
                              PYTHON_FILES ${PYTHON_FILES})
endif()
if (CARIBOU_WITH_PYTHON_3)
    caribou_add_python_module(CaribouPython
                              PYTHON_VERSION 3
                              DESTINATION Caribou
                              PYTHON_FILES ${PYTHON_FILES})
endif()