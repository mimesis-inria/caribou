set(DIRECTORY Caribou)

set(PYTHON_FILES
    __init__.py
)

set(PYTHON_TEST_FILES
)

set(HEADER_FILES
    Caribou.h
)

set(SOURCE_FILES
    Caribou.cpp
)

if (CARIBOU_WITH_PYTHON_3)
    caribou_add_python_module(CaribouPython
            DESTINATION ${DIRECTORY}
            SOURCE_FILES ${SOURCE_FILES}
            PYTHON_FILES ${PYTHON_FILES}
            PYTHON_TEST_FILES ${PYTHON_TEST_FILES}
    )
endif()