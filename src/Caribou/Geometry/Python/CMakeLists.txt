set(DIRECTORY Caribou/Geometry)

set(PYTHON_FILES
    __init__.py
)

set(PYTHON_TEST_FILES
        test/Caribou.Geometry.Test.py
)

set(HEADER_FILES
        Element.h
)

set(SOURCE_FILES
        Geometry.cpp
        Hexahedron.cpp
        Quad.cpp
        RectangularHexahedron.cpp
        Segment.cpp
        Tetrahedron.cpp
        Triangle.cpp
)

unset(PYTHON_MODULE_PREFIX)
unset(PYTHON_MODULE_EXTENSION)

if (CARIBOU_WITH_PYTHON_3)
    caribou_add_python_module(CaribouGeometryPython
                              DESTINATION ${DIRECTORY}
                              SOURCE_FILES ${SOURCE_FILES}
                              PYTHON_FILES ${PYTHON_FILES}
                              PYTHON_TEST_FILES ${PYTHON_TEST_FILES}
                              DEPENDS Geometry)
endif()
