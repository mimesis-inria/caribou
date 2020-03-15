# Configure and install Caribou's main header files
set(header_files
    config.h.in
    macros.h
    Traits.h
)

foreach(header_file ${header_files})
    string(REGEX REPLACE "\\.in$" "" filename ${header_file})
    configure_file("${header_file}" "${CMAKE_CURRENT_BINARY_DIR}/${filename}")
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${filename} DESTINATION "include/Caribou" COMPONENT headers)
endforeach()

# Caribou's component
set(Caribou_COMPONENTS 
    Algebra
    Geometry
    Mechanics
    Topology
)


foreach(component ${Caribou_COMPONENTS})
    add_subdirectory(${component})
endforeach()

configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/CaribouConfig.cmake.in"
    "${CMAKE_BINARY_DIR}/cmake/Caribou/CaribouConfig.cmake"
    INSTALL_DESTINATION
    lib/cmake/Caribou
)
write_basic_package_version_file(
        ${CMAKE_BINARY_DIR}/cmake/Caribou/CaribouConfigVersion.cmake
        VERSION ${Caribou_VERSION}
        COMPATIBILITY AnyNewerVersion)

install(
    EXPORT CaribouTargets 
    NAMESPACE Caribou::
    DESTINATION "lib/cmake/Caribou" 
    COMPONENT headers)
install(FILES 
        "${CMAKE_BINARY_DIR}/cmake/Caribou/CaribouConfig.cmake"
        "${CMAKE_BINARY_DIR}/cmake/Caribou/CaribouConfigVersion.cmake"
        DESTINATION 
        "lib/cmake/Caribou")
