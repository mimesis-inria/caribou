project(Caribou.unittests.Topology)

set(HEADER_FILES
    Grid/Grid_1D.h
    Grid/Grid_2D.h
    Grid/Grid_3D.h
    topology_test.h
)

set(SOURCE_FILES
    Grid/Grid.cpp
    test_barycentric_container.cpp
    test_domain.cpp
    test_mesh.cpp
    main.cpp
)

set(MESH_FILES
    meshes/1D_linear.vtk
    meshes/1D_quadratic.vtk
    meshes/2D_quad_linear.vtk
    meshes/2D_quad_quadratic.vtk
    meshes/2D_triangle_linear.vtk
    meshes/2D_triangle_quadratic.vtk
    meshes/3D_hexahedron_linear.vtk
    meshes/3D_hexahedron_quadratic.vtk
    meshes/3D_quad_linear.vtk
    meshes/3D_quad_quadratic.vtk
    meshes/3D_tetrahedron_linear.vtk
    meshes/3D_tetrahedron_quadratic.vtk
    meshes/3D_triangle_linear.vtk
    meshes/3D_triangle_quadratic.vtk
    meshes/deformed_liver_surface.stl
    meshes/deformed_liver_volume_hexahedrons.vtu
    meshes/deformed_liver_volume_tetrahedrons.vtu
)

if (CARIBOU_WITH_VTK)
    list(APPEND SOURCE_FILES test_vtkreader.cpp)
endif()

enable_testing()

add_executable(${PROJECT_NAME} ${SOURCE_FILES})
target_link_libraries(${PROJECT_NAME} PUBLIC Caribou::Topology)
target_link_libraries(${PROJECT_NAME} PUBLIC gtest)

if (NOT WIN32)
    target_link_libraries(${PROJECT_NAME} PUBLIC pthread)
endif()


list(APPEND target_rpath
    "$ORIGIN/../lib"
    "$ORIGIN/../../../lib"
    "@executable_path/../lib"
    "@executable_path/../../../lib"
)

set_target_properties(${PROJECT_NAME} PROPERTIES INSTALL_RPATH "${target_rpath}" )

target_link_libraries(${PROJECT_NAME} PUBLIC
                      $<$<AND:$<PLATFORM_ID:Linux>,$<VERSION_LESS:$<CXX_COMPILER_VERSION>,8.0>>:stdc++fs>
                      $<$<AND:$<PLATFORM_ID:Darwin>,$<VERSION_LESS:$<CXX_COMPILER_VERSION>,9.0>>:c++fs>)
target_compile_definitions(${PROJECT_NAME} PUBLIC
                           $<$<AND:$<PLATFORM_ID:Linux>,$<VERSION_LESS:$<CXX_COMPILER_VERSION>,8.0>>:LEGACY_CXX>
                           $<$<AND:$<PLATFORM_ID:Darwin>,$<VERSION_LESS:$<CXX_COMPILER_VERSION>,9.0>>:LEGACY_CXX>)

foreach(file ${MESH_FILES})
    file(RELATIVE_PATH path_from_package "${CMAKE_CURRENT_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}/${file}")
    get_filename_component(dir_from_package ${path_from_package} DIRECTORY)
    configure_file("${file}" "${CMAKE_BINARY_DIR}/bin/${path_from_package}" COPYONLY)
    install(FILES "${CMAKE_BINARY_DIR}/bin/${path_from_package}" DESTINATION "bin/${dir_from_package}")
endforeach()

install(
    TARGETS ${PROJECT_NAME}
    EXPORT Caribou
)
