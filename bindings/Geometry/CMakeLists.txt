project(Bindings.Geometry)

set(PYTHON_FILES
    __init__.py
)

set(PYTHON_TEST_FILES
    pytest/Caribou.Geometry.Test.py
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
                          TARGET_ALIAS      Caribou::${PROJECT_NAME}
                          TARGET_DEPENDS    Caribou::Geometry Caribou::Bindings.Config
                          COMPONENT_NAME    Bindings
                          PACKAGE_NAME      Caribou
                          DESTINATION       Caribou/Geometry
                          SOURCE_FILES      ${SOURCE_FILES}
                          PUBLIC_HEADERS    ${HEADER_FILES}
                          HEADER_PREFIX     Caribou/Bindings/Geometry
                          PYTHON_FILES      ${PYTHON_FILES}
                          PYTHON_TEST_FILES ${PYTHON_TEST_FILES}
)
