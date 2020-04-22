project(Caribou.Topology.Test)

set(HEADER_FILES
        Grid/Grid_1D.h
        Grid/Grid_2D.h
        Grid/Grid_3D.h
        topology_test.h
)

set(SOURCE_FILES
        Grid/Grid.cpp
        test_mesh.cpp
        test_domain.cpp
        main.cpp)

set(MESH_FILES
        meshes/1D_linear.vtk
        meshes/1D_quadratic.vtk
        meshes/2D_quad_linear.vtk
        meshes/2D_quad_quadratic.vtk
        meshes/2D_triangle_linear.vtk
        meshes/2D_triangle_quadratic.vtk
        meshes/3D_quad_linear.vtk
        meshes/3D_quad_quadratic.vtk
        meshes/3D_tetrahedron_linear.vtk
        meshes/3D_tetrahedron_quadratic.vtk
        meshes/3D_triangle_linear.vtk
        meshes/3D_triangle_quadratic.vtk
)

if (CARIBOU_WITH_VTK)
    list(APPEND SOURCE_FILES test_vtkreader.cpp)
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