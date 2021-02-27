project(Topology)

set(HEADER_FILES
    config.h.in
    BarycentricContainer.h
    BaseMesh.h
    BaseDomain.h
    Domain.h
    Grid/Grid.h
    Grid/Internal/BaseGrid.h
    Grid/Internal/BaseMultidimensionalGrid.h
    Grid/Internal/BaseUnidimensionalGrid.h
    HashGrid.h
    Mesh.h
)

set(TARGET_TYPE "INTERFACE")
set(TARGET_VISIBILITY "INTERFACE")
if (CARIBOU_WITH_VTK)
    find_package(VTK COMPONENTS ${CARIBOU_VTK_MODULES} REQUIRED)
    set(TARGET_TYPE "SHARED")
    set(TARGET_VISIBILITY "PUBLIC")
    list(APPEND HEADER_FILES IO/VTKReader.h)
    list(APPEND SOURCE_FILES IO/VTKReader.cpp)
    if (VTK_VERSION VERSION_LESS "8.90.0")
        # old system
        include(${VTK_USE_FILE})
    endif()
endif()

find_package(Eigen3 QUIET REQUIRED)

add_library(${PROJECT_NAME} ${TARGET_TYPE} ${SOURCE_FILES})
add_library(Caribou::${PROJECT_NAME} ALIAS ${PROJECT_NAME})

target_link_libraries (${PROJECT_NAME} ${TARGET_VISIBILITY} Eigen3::Eigen)
target_link_libraries(${PROJECT_NAME} ${TARGET_VISIBILITY} Caribou::Geometry)
target_link_libraries(${PROJECT_NAME} ${TARGET_VISIBILITY} $<$<AND:$<CXX_COMPILER_ID:GNU>,$<VERSION_LESS:$<CXX_COMPILER_VERSION>,9.0>>:stdc++fs>)
if (CARIBOU_WITH_VTK)
    target_link_libraries(${PROJECT_NAME} PUBLIC ${VTK_LIBRARIES})
endif()
if (VTK_VERSION VERSION_GREATER_EQUAL "8.90.0")
    vtk_module_autoinit(
        TARGETS ${PROJECT_NAME}
        MODULES ${VTK_LIBRARIES}
    )
endif()

# Add the target to the component of the same name. This will enable to do:
# CMakeLists.txt
#    find_package(Caribou COMPONENTS Topology)
#    add_library(my_cool_library SHARED)
#    target_link_libraries(my_cool_library PRIVATE Caribou::Topology)
caribou_add_target_to_component (
    TARGET_NAME    ${PROJECT_NAME}
    COMPONENT_NAME ${PROJECT_NAME}
    PACKAGE_NAME   Caribou
    PUBLIC_HEADERS ${HEADER_FILES}
)

# This is required to be able to do "find_package(Caribou COMPONENTS Topology)"
caribou_add_component_to_package (
    COMPONENT_NAME ${PROJECT_NAME}
    PACKAGE_NAME   Caribou
)