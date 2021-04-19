add_subdirectory(Config)
add_subdirectory(Geometry)
add_subdirectory(Topology)

caribou_add_component_to_package(
    COMPONENT_NAME Python
    PACKAGE_NAME   Caribou
)