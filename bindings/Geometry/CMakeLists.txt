project(Caribou.Bindings.Geometry)

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
                          TARGET_ALIAS      Caribou::Bindings.Geometry
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

if (NOT WIN32)
    # Avoid an annoying warning since the bindings code are in "caribou/bindings/*" but the headers are installed into
    # "Caribou/Bindings/*" to be consistent with SOFA and the rest of caribou libraries.
    target_compile_options(${PROJECT_NAME} PRIVATE "-Wno-nonportable-include-path")
endif()