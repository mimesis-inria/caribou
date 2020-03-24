project(Geometry)

set(HEADER_FILES
    Hexahedron.h
    Quad.h
    RectangularQuad.h
    RectangularHexahedron.h
    Segment.h
    Traits.h
    Triangle.h
    Tetrahedron.h
    Internal/BaseHexahedron.h
    Internal/BaseQuad.h
    Internal/BaseSegment.h
    Internal/BaseTetrahedron.h
    Internal/BaseTriangle.h
    Interpolation/CanonicalElement.h
    Interpolation/Hexahedron.h
    Interpolation/Quad.h
    Interpolation/Tetrahedron.h
    Interpolation/Triangle.h
    Interpolation/Segment.h)


find_package(Eigen3 REQUIRED NO_MODULE)

add_library(${PROJECT_NAME} INTERFACE)
add_library(Caribou::${PROJECT_NAME} ALIAS ${PROJECT_NAME})
target_include_directories(${PROJECT_NAME} INTERFACE "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../../>")
target_include_directories(${PROJECT_NAME} INTERFACE "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/../../>")
target_include_directories(${PROJECT_NAME} INTERFACE "$<INSTALL_INTERFACE:include>")
target_link_libraries (${PROJECT_NAME} INTERFACE Eigen3::Eigen)
target_compile_features(${PROJECT_NAME} INTERFACE cxx_std_14)

caribou_install_target(Caribou ${PROJECT_NAME} ${HEADER_FILES})

if (CARIBOU_BUILD_TESTS)
    add_subdirectory(test)
endif()

if (CARIBOU_WITH_PYTHON_2 OR CARIBOU_WITH_PYTHON_3)
    add_subdirectory(Python)
endif()
