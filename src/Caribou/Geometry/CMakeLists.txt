project(Geometry)

set(HEADER_FILES
    Base/BaseHexahedron.h
    Base/BaseQuad.h
    Base/BaseRectangularHexahedron.h
    Base/BaseRectangularQuad.h
    Base/BaseSegment.h
    Base/BaseTetrahedron.h
    Base/BaseTriangle.h
    Element.h
    Hexahedron.h
    Hexahedron_FEniCS.h
    Quad.h
    RectangularHexahedron.h
    RectangularQuad.h
    Segment.h
    Tetrahedron.h
    Triangle.h
)

find_package(Eigen3 QUIET REQUIRED)

add_library(${PROJECT_NAME} INTERFACE)
add_library(Caribou::${PROJECT_NAME} ALIAS ${PROJECT_NAME})

target_link_libraries (${PROJECT_NAME} INTERFACE Caribou::Config)
target_link_libraries (${PROJECT_NAME} INTERFACE Eigen3::Eigen)

# Add the target Geometry to the component of the same name. This will enable to do:
# CMakeLists.txt
#    find_package(Caribou COMPONENTS Geometry)
#    add_library(my_cool_library SHARED)
#    target_link_libraries(my_cool_library PRIVATE Caribou::Geometry)
caribou_add_target_to_component (
    TARGET_NAME    ${PROJECT_NAME}
    COMPONENT_NAME ${PROJECT_NAME}
    PACKAGE_NAME   Caribou
    PUBLIC_HEADERS ${HEADER_FILES}
)

# This is required to be able to do "find_package(Caribou COMPONENTS Geometry)"
caribou_add_component_to_package (
    COMPONENT_NAME ${PROJECT_NAME}
    PACKAGE_NAME   Caribou
)
