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

add_library(${PROJECT_NAME} ${TARGET_TYPE} ${SOURCE_FILES})
add_library(Caribou::${PROJECT_NAME} ALIAS ${PROJECT_NAME})

target_include_directories(${PROJECT_NAME} ${TARGET_VISIBILITY} "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../../>")
target_include_directories(${PROJECT_NAME} ${TARGET_VISIBILITY} "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/../../>")
target_include_directories(${PROJECT_NAME} ${TARGET_VISIBILITY} "$<INSTALL_INTERFACE:include>")
target_link_libraries(${PROJECT_NAME} ${TARGET_VISIBILITY} Caribou::Geometry)
target_compile_features(${PROJECT_NAME} ${TARGET_VISIBILITY} cxx_std_17)
target_link_libraries(${PROJECT_NAME} ${TARGET_VISIBILITY} $<$<AND:$<CXX_COMPILER_ID:GNU>,$<VERSION_LESS:$<CXX_COMPILER_VERSION>,9.0>>:stdc++fs>)

caribou_install_target(Caribou ${PROJECT_NAME} ${HEADER_FILES})

if (CARIBOU_WITH_VTK)
    target_link_libraries(${PROJECT_NAME} PUBLIC ${VTK_LIBRARIES})
endif()

if (CARIBOU_BUILD_TESTS)
    add_subdirectory(test)
endif()

if (CARIBOU_WITH_PYTHON_3)
    add_subdirectory(Python)
endif()
