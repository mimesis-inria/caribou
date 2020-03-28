project(Geometry)

set(HEADER_FILES
        BaseHexahedron.h
        BaseSegment.h
        BaseQuad.h
        BaseRectangularQuad.h
        BaseTetrahedron.h
        BaseTriangle.h
        Element.h
        RectangularQuad.h
        Segment.h
        Quad.h
        Tetrahedron.h
        Triangle.h)


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

#if (CARIBOU_WITH_PYTHON_2 OR CARIBOU_WITH_PYTHON_3)
#    add_subdirectory(Python)
#endif()