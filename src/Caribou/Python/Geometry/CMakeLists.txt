project(Caribou.Python.Geometry)

set(PYTHON_FILES
    __init__.py
)

set(PYTHON_TEST_FILES
    pytest/Caribou_Geometry_Test.py
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

caribou_add_python_module(Geometry
                          TARGET_NAME       ${PROJECT_NAME}
                          TARGET_ALIAS      Caribou::Python.Geometry
                          TARGET_DEPENDS    Caribou::Geometry Caribou::Python.Config
                          COMPONENT_NAME    Python
                          PACKAGE_NAME      Caribou
                          DESTINATION       Caribou/Geometry
                          SOURCE_FILES      ${SOURCE_FILES}
                          PUBLIC_HEADERS    ${HEADER_FILES}
                          HEADER_SRC_PREFIX "${CMAKE_CURRENT_SOURCE_DIR}/../../.."
                          PYTHON_FILES      ${PYTHON_FILES}
                          PYTHON_TEST_FILES ${PYTHON_TEST_FILES}
)

