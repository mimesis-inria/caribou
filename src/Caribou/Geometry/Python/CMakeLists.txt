set(DIRECTORY Caribou/Geometry)

set(PYTHON_FILES
    __init__.py
)

set(PYTHON_TEST_FILES
    tests/test_segment.py
)

set(HEADER_FILES
    Hexahedron.h
    Segment.h
    Tetrahedron.h
)

set(SOURCE_FILES
    Geometry.cpp
    Hexahedron.cpp
    Segment.cpp
    Tetrahedron.cpp
)

unset(PYTHON_MODULE_PREFIX)
unset(PYTHON_MODULE_EXTENSION)

if (CARIBOU_WITH_PYTHON_2)
    caribou_add_python_module(CaribouGeometryPython
                              PYTHON_VERSION 2.7
                              DESTINATION Caribou/Geometry
                              SOURCE_FILES ${SOURCE_FILES}
                              PYTHON_FILES ${PYTHON_FILES}
                              PYTHON_TEST_FILES ${PYTHON_TEST_FILES}
                              DEPENDS Geometry)
endif()

if (CARIBOU_WITH_PYTHON_3)
    caribou_add_python_module(CaribouGeometryPython
                              PYTHON_VERSION 3
                              DESTINATION Caribou/Geometry
                              SOURCE_FILES ${SOURCE_FILES}
                              PYTHON_FILES ${PYTHON_FILES}
                              PYTHON_TEST_FILES ${PYTHON_TEST_FILES}
                              DEPENDS Geometry)
endif()
