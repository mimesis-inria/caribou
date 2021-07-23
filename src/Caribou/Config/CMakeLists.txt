project(Config)

set (HEADER_FILES
     config.h.in
     constants.h.in
     macros.h.in
     traits.h.in
     version.h.in
)

add_library(${PROJECT_NAME} INTERFACE)
add_library(Caribou::${PROJECT_NAME} ALIAS ${PROJECT_NAME})

set(CARIBOU_VERSION_STR "\"${CARIBOU_VERSION}\"")
string(REPLACE "." "" CARIBOU_VERSION_NUM ${CARIBOU_VERSION})
string(REPLACE " " "" CARIBOU_VERSION_NUM ${CARIBOU_VERSION_NUM})

# Add the target Config to the component of the same name. This will enable to do:
# CMakeLists.txt
#    find_package(Caribou COMPONENTS Config)
#    add_library(my_cool_library SHARED)
#    target_link_libraries(my_cool_library PRIVATE Caribou::Config)
caribou_add_target_to_component (
    TARGET_NAME    ${PROJECT_NAME}
    COMPONENT_NAME ${PROJECT_NAME}
    PACKAGE_NAME   Caribou
    PUBLIC_HEADERS ${HEADER_FILES}
    HEADER_BUILD_PREFIX  "${CMAKE_BINARY_DIR}/include/Caribou"
    HEADER_INSTALL_PREFIX  "include/Caribou"
)

# This is required to be able to do "find_package(Caribou COMPONENTS Config)"
caribou_add_component_to_package (
    COMPONENT_NAME ${PROJECT_NAME}
    PACKAGE_NAME   Caribou
)