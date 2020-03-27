project(Caribou.Topology.Test)

set(HEADER_FILES
        Grid/Grid.h
        Grid/Grid_1D.h
        Grid/Grid_2D.h
        Grid/Grid_3D.h
        test_unstructured_grid.h
        test_vtkreader.h
)

set(SOURCE_FILES
        main.cpp)

set(MESH_FILES
        meshes/1D_linear.vtk)

if (CARIBOU_WITH_VTK)
    list(APPEND ${HEADER_FILES} test_vtkreader.h)
endif()

enable_testing()
find_package(GTest QUIET REQUIRED)

add_executable(${PROJECT_NAME} ${SOURCE_FILES})
target_link_libraries(${PROJECT_NAME} ${GTEST_BOTH_LIBRARIES} Topology)
target_include_directories(${PROJECT_NAME} PUBLIC "$<INSTALL_INTERFACE:include>")
target_include_directories(${PROJECT_NAME} PUBLIC "$<BUILD_INTERFACE:${GTEST_INCLUDE_DIR}>")

foreach(file ${MESH_FILES})
    file(RELATIVE_PATH path_from_package "${CMAKE_CURRENT_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}/${file}")
    configure_file("${file}" "${CMAKE_BINARY_DIR}/bin/${path_from_package}" COPYONLY)
endforeach()