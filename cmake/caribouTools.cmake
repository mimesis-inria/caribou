function(dump_cmake_variables)
    get_cmake_property(_variableNames VARIABLES)
    list (SORT _variableNames)
    foreach (_variableName ${_variableNames})
        if (ARGV0)
            unset(MATCHED)
            string(REGEX MATCH ${ARGV0} MATCHED ${_variableName})
            if (NOT MATCHED)
                continue()
            endif()
        endif()
        message(STATUS "${_variableName}=${${_variableName}}")
    endforeach()
endfunction()

function(caribou_link from to)
    IF(WIN32)
        SET(LINK copy_if_different)
    ELSE(WIN32)
        SET(LINK create_symlink)
    ENDIF(WIN32)

    execute_process(COMMAND ${CMAKE_COMMAND} -E ${LINK} ${from} ${to})
endfunction()

# caribou_add_python_module ( <name>
#                             TARGET_NAME          <target_name>        # Default to <name>
#                             TARGET_ALIAS         <target_alias>
#                             TARGET_VERSION       <version>            # defaults to ${PACKAGE_NAME}_VERSION
#                             TARGET_DEPENDS       [<target1, ...>]
#                             COMPONENT_NAME       [<component_name>]
#                             PACKAGE_NAME         [<package_name>]
#                             DESTINATION          [<package_name>]
#                             PREFIX               [<prefix>]           # Default to python3/site-packages
#                             SOURCE_FILES         [<file1, ...>]
#                             PUBLIC_HEADERS       [<file1, ...>]
#                             PRIVATE_HEADERS      [<file1, ...>]
#                             HEADER_SRC_PREFIX     <[relative_path]>   # Path from where the relative path to the header files will be
#                                                                       # computed. Default to "${CMAKE_CURRENT_SOURCE_DIR}/../..".
#                             HEADER_BUILD_PREFIX   <[relative_path]>   # Path where the header file will be configured. Default to "${CMAKE_BINARY_DIR}/include/__P__"
#                                                                       # where __P__ is the relative path between ${HEADER_SRC_PREFIX} and the header file path.
#                             HEADER_INSTALL_PREFIX <[relative_path]>   # Path where the header file will be installed. Default to "include/__P__".
#                                                                       # where __P__ is the relative path between ${HEADER_SRC_PREFIX} and the header file path.
#                             PYTHON_FILES         [<file1, ...>]
#                             PYTHON_TEST_FILES    [<file1, ...>]
#                             LIBRARY_DESTINATION  [<folder_path>]      # default to "lib/${PREFIX}/${DESTINATION}"
#                             RUNTIME_DESTINATION  [<folder_path>]      # default to "bin"
# )
#
# Creates a python modules from ${SOURCE_FILES} and deploy the resulting binary named "${NAME}.so" into
# the folder ${LIBRARY_DESTINATION}. The ${PUBLIC_HEADERS} can be used to specify which header files
# should be installed. Both ${PUBLIC_HEADERS} and ${PRIVATE_HEADERS} will be
# added to the compilation header search path. In addition, ${PYTHON_FILES} will be configured (if
# suffixed by .in) or linked into the compilation directory ${CMAKE_BINARY_DIR}/${LIBRARY_DESTINATION}
# and then installed in ${CMAKE_INSTALL_PREFIX}/${LIBRARY_DESTINATION}. Finally, ${PYTHON_TEST_FILES} will
# be configured (if suffixed by .in) or linked into the compilation directory ${CMAKE_BINARY_DIR}/${RUNTIME_DESTINATION}
# and then installed in ${CMAKE_INSTALL_PREFIX}/${RUNTIME_DESTINATION}.
function(caribou_add_python_module NAME)
    set(options QUIET)
    set(oneValueArgs TARGET_NAME TARGET_ALIAS TARGET_VERSION COMPONENT_NAME PACKAGE_NAME DESTINATION PREFIX HEADER_SRC_PREFIX HEADER_BUILD_PREFIX HEADER_INSTALL_PREFIX LIBRARY_DESTINATION RUNTIME_DESTINATION)
    set(multiValueArgs SOURCE_FILES PUBLIC_HEADERS PRIVATE_HEADERS PYTHON_FILES PYTHON_TEST_FILES TARGET_DEPENDS)

    cmake_parse_arguments(A "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

    set(PYBIND11_PYTHON_VERSION 3)
    set(PYBIND11_FINDPYTHON ON)
    find_package(Python ${PYBIND11_PYTHON_VERSION} COMPONENTS Interpreter Development QUIET)
    find_package(pybind11 QUIET CONFIG REQUIRED)

    set(TARGET_NAME     ${A_TARGET_NAME})
    set(TARGET_ALIAS    ${A_TARGET_ALIAS})
    set(TARGET_VERSION  ${A_TARGET_VERSION})
    set(TARGET_DEPENDS  ${A_TARGET_DEPENDS})
    set(COMPONENT_NAME  ${A_COMPONENT_NAME})
    set(PACKAGE_NAME    ${A_PACKAGE_NAME})
    set(DESTINATION     ${A_DESTINATION})
    set(PREFIX          ${A_PREFIX})
    set(SOURCE_FILES      "${A_SOURCE_FILES}")
    set(PUBLIC_HEADERS    "${A_PUBLIC_HEADERS}")
    set(PRIVATE_HEADERS   "${A_PRIVATE_HEADERS}")
    set(HEADER_SRC_PREFIX     "${A_HEADER_SRC_PREFIX}")
    set(HEADER_BUILD_PREFIX   "${A_HEADER_BUILD_PREFIX}")
    set(HEADER_INSTALL_PREFIX "${A_HEADER_INSTALL_PREFIX}")
    set(PYTHON_FILES      "${A_PYTHON_FILES}")
    set(PYTHON_TEST_FILES "${A_PYTHON_TEST_FILES}")
    set(LIBRARY_DESTINATION ${A_LIBRARY_DESTINATION})
    set(RUNTIME_DESTINATION ${A_RUNTIME_DESTINATION})

    if (NOT TARGET_NAME)
        set(TARGET_NAME "${NAME}")
    endif()

    if (NOT PREFIX)
        set(PREFIX "python3/site-packages")
    endif()

    if (NOT TARGET_VERSION)
        set(TARGET_VERSION ${${PACKAGE_NAME}_VERSION})
    endif()

    if (NOT LIBRARY_DESTINATION)
        set(LIBRARY_DESTINATION "lib/${PREFIX}/${DESTINATION}")
    endif()

    if (NOT RUNTIME_DESTINATION)
        set(RUNTIME_DESTINATION "bin")
    endif()

    if (A_SOURCE_FILES)
        set(PYBIND11_CPP_STANDARD -std=c++17)

        # We are doing manually what's usually done with pybind11_add_module(${TARGET_NAME} SHARED "${A_SOURCE_FILES}")
        # since we got some problems on MacOS using recent versions of pybind11 where the SHARED argument wasn't taken
        # into account
        python_add_library(${TARGET_NAME} SHARED "${A_SOURCE_FILES}")

        if (A_TARGET_ALIAS)
            add_library(${A_TARGET_ALIAS} ALIAS ${TARGET_NAME})
        endif ()

        if ("${pybind11_VERSION}" VERSION_GREATER_EQUAL "2.6.0")
            target_link_libraries(${TARGET_NAME} PUBLIC pybind11::headers)
            target_link_libraries(${TARGET_NAME} PRIVATE pybind11::embed)
            target_link_libraries(${TARGET_NAME} PRIVATE pybind11::lto)
            if(MSVC)
                target_link_libraries(${TARGET_NAME} PRIVATE pybind11::windows_extras)
            endif()

            pybind11_extension(${TARGET_NAME})
            pybind11_strip(${TARGET_NAME})
        else()
            target_link_libraries(${TARGET_NAME} PUBLIC pybind11::module)

            # Equivalent to pybind11_extension(${TARGET_NAME}) which doesn't exists on pybind11 versions < 5
            set_target_properties(${TARGET_NAME} PROPERTIES PREFIX "" SUFFIX "${PYTHON_MODULE_EXTENSION}")

            if(NOT MSVC AND NOT ${CMAKE_BUILD_TYPE} MATCHES Debug|RelWithDebInfo)
                # Equivalent to pybind11_strip(${TARGET_NAME}) which doesn't exists on pybind11 versions < 5
                # Strip unnecessary sections of the binary on Linux/macOS
                if(CMAKE_STRIP)
                    if(APPLE)
                        set(x_opt -x)
                    endif()

                    add_custom_command(
                        TARGET ${TARGET_NAME}
                        POST_BUILD
                        COMMAND ${CMAKE_STRIP} ${x_opt} $<TARGET_FILE:${TARGET_NAME}>)
                endif()
            endif()
        endif()

        target_link_libraries(${TARGET_NAME} PUBLIC ${TARGET_DEPENDS})

        if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
            target_compile_options(${TARGET_NAME} PRIVATE -fsized-deallocation)
        endif()

        target_compile_features(${TARGET_NAME} PUBLIC cxx_std_17)

        set_target_properties(
            ${TARGET_NAME}
            PROPERTIES
                OUTPUT_NAME ${NAME}
                LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${LIBRARY_DESTINATION}"
        )
        message(STATUS "${TARGET_NAME}: ${CMAKE_BINARY_DIR}/${LIBRARY_DESTINATION}")

        caribou_add_target_to_component(
            TARGET_NAME            ${TARGET_NAME}
            COMPONENT_NAME         ${COMPONENT_NAME}
            PACKAGE_NAME           ${PACKAGE_NAME}
            TARGET_VERSION         ${TARGET_VERSION}
            PUBLIC_HEADERS         ${PUBLIC_HEADERS}
            PRIVATE_HEADERS        ${PRIVATE_HEADERS}
            HEADER_SRC_PREFIX      ${HEADER_SRC_PREFIX}
            HEADER_BUILD_PREFIX    ${HEADER_BUILD_PREFIX}
            HEADER_INSTALL_PREFIX  ${HEADER_INSTALL_PREFIX}
            LIBRARY_DESTINATION    ${LIBRARY_DESTINATION}
            RUNTIME_DESTINATION    ${RUNTIME_DESTINATION}
        )

        # This will get the relative path from the current binding library to the plugin's "lib" directory.
        # As an example, build/lib/python3.7/site-packages/Caribou/Geometry/CaribouGeometryPython.cpython-39-x86_64-linux.so
        # will give "../../../.."
        # We can thereby add this path to the installation RPATH in order to let the binding library find Caribou's
        # libraries.
        file(RELATIVE_PATH path_to_lib  "${CMAKE_BINARY_DIR}/${LIBRARY_DESTINATION}" "${CMAKE_BINARY_DIR}/lib")
        list(APPEND target_rpath
            "$ORIGIN/${path_to_lib}"
            "$loader_path/${path_to_lib}"
        )
        set_target_properties(
            ${TARGET_NAME}
            PROPERTIES
                INSTALL_RPATH "${target_rpath}"
        )
    endif ()

    foreach(t ${PYTHON_FILES})
        file(RELATIVE_PATH path_from_current "${CMAKE_CURRENT_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}/${t}")
        get_filename_component(dir_from_current ${path_from_current} DIRECTORY)
        get_filename_component(output_filename ${path_from_current} NAME)

        file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/${LIBRARY_DESTINATION}/${dir_from_current}")

        if (output_filename MATCHES "\\.in$")
            string(REGEX REPLACE "\\.in$" "" output_filename ${output_filename})
            configure_file("${t}" "${CMAKE_BINARY_DIR}/${LIBRARY_DESTINATION}/${dir_from_current}/${output_filename}")
        else()
            caribou_link("${CMAKE_CURRENT_SOURCE_DIR}/${dir_from_current}/${output_filename}" "${CMAKE_BINARY_DIR}/${LIBRARY_DESTINATION}/${dir_from_current}/${output_filename}")
        endif()

        install(FILES "${CMAKE_BINARY_DIR}/${LIBRARY_DESTINATION}/${dir_from_current}/${output_filename}" DESTINATION ${LIBRARY_DESTINATION}/${dir_from_current})
    endforeach()

    foreach(t ${PYTHON_TEST_FILES})
        file(RELATIVE_PATH path_from_current "${CMAKE_CURRENT_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}/${t}")
        get_filename_component(dir_from_current ${path_from_current} DIRECTORY)
        get_filename_component(output_filename ${path_from_current} NAME)

        file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/${RUNTIME_DESTINATION}/${dir_from_current}")

        if (output_filename MATCHES "\\.in$")
            string(REGEX REPLACE "\\.in$" "" output_filename ${output_filename})
            configure_file("${t}" "${CMAKE_BINARY_DIR}/${RUNTIME_DESTINATION}/${dir_from_current}/${output_filename}")
        else()
            caribou_link("${CMAKE_CURRENT_SOURCE_DIR}/${dir_from_current}/${output_filename}" "${CMAKE_BINARY_DIR}/${RUNTIME_DESTINATION}/${dir_from_current}/${output_filename}")
        endif()

        install(FILES "${CMAKE_BINARY_DIR}/${RUNTIME_DESTINATION}/${dir_from_current}/${output_filename}" DESTINATION ${RUNTIME_DESTINATION}/${dir_from_current})
    endforeach()

endfunction()

# caribou_create_package (
#     PACKAGE_NAME <package_name>
# )
#
# Create a CMake package that will be later found by find_package(PACKAGE_NAME). This macro handles the creation
# and installation of PACKAGE_NAMEVersion.cmake and PACKAGE_NAMEConfig.cmake (if found).
macro(caribou_create_package PACKAGE_NAME PACKAGE_VERSION)

    if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${PACKAGE_NAME}Config.cmake.in")
        get_property(component_names GLOBAL PROPERTY "${PACKAGE_NAME}_COMPONENTS")
        set(${PACKAGE_NAME}_COMPONENTS ${component_names})
        configure_package_config_file (
            "${CMAKE_CURRENT_SOURCE_DIR}/${PACKAGE_NAME}Config.cmake.in"
            "${CMAKE_BINARY_DIR}/cmake/${PACKAGE_NAME}/${PACKAGE_NAME}Config.cmake"
            INSTALL_DESTINATION
                lib/cmake/${PACKAGE_NAME}
        )

        write_basic_package_version_file (
            "${CMAKE_BINARY_DIR}/cmake/${PACKAGE_NAME}/${PACKAGE_NAME}ConfigVersion.cmake"
            VERSION ${PACKAGE_VERSION}
            COMPATIBILITY ExactVersion
        )

        install (
            FILES
                "${CMAKE_BINARY_DIR}/cmake/${PACKAGE_NAME}/${PACKAGE_NAME}Config.cmake"
                "${CMAKE_BINARY_DIR}/cmake/${PACKAGE_NAME}/${PACKAGE_NAME}ConfigVersion.cmake"
            DESTINATION
                "lib/cmake/${PACKAGE_NAME}"
        )
    else()
        message(FATAL_ERROR "${CMAKE_CURRENT_SOURCE_DIR}/${PACKAGE_NAME}Config.cmake.in does not exist")
    endif()

endmacro()

# caribou_add_component_to_package (
#     PACKAGE_NAME       <package_name>
#     COMPONENT_NAME     <component_name>
#     COMPONENT_VERSION  <component_version>  # Default to ${PACKAGE_NAME}_VERSION
# )
#
# Create a CMake component and add it to a given package. This component will be later found using
# find_package(PACKAGE_NAME). This macro handles the creation and installation of PACKAGE_NAMEVersion.cmake
# and PACKAGE_NAMEConfig.cmake (if found).
macro(caribou_add_component_to_package)
    set(oneValueArgs COMPONENT_NAME PACKAGE_NAME)
    set(multiValueArgs )
    set(optionalArgs )
    cmake_parse_arguments("ARG" "${optionalArgs}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

    set(COMPONENT_NAME ${ARG_COMPONENT_NAME})
    set(PACKAGE_NAME   ${ARG_PACKAGE_NAME})

    # Target version
    set(COMPONENT_VERSION ${ARG_COMPONENT_VERSION})
    if (NOT COMPONENT_VERSION)
        set(COMPONENT_VERSION ${${PACKAGE_NAME}_VERSION})
    endif()

    if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${COMPONENT_NAME}Config.cmake.in")
        configure_package_config_file (
            "${CMAKE_CURRENT_SOURCE_DIR}/${COMPONENT_NAME}Config.cmake.in"
            "${CMAKE_BINARY_DIR}/cmake/${PACKAGE_NAME}/${COMPONENT_NAME}Config.cmake"
            INSTALL_DESTINATION
                lib/cmake/${PACKAGE_NAME}
        )

        write_basic_package_version_file (
            "${CMAKE_BINARY_DIR}/cmake/${PACKAGE_NAME}/${COMPONENT_NAME}ConfigVersion.cmake"
            VERSION ${COMPONENT_VERSION}
            COMPATIBILITY ExactVersion
        )

        install (
            FILES
            "${CMAKE_BINARY_DIR}/cmake/${PACKAGE_NAME}/${COMPONENT_NAME}Config.cmake"
            "${CMAKE_BINARY_DIR}/cmake/${PACKAGE_NAME}/${COMPONENT_NAME}ConfigVersion.cmake"
            DESTINATION
            "lib/cmake/${PACKAGE_NAME}"
        )

        install (
            EXPORT ${COMPONENT_NAME}Targets
            DESTINATION "lib/cmake/${PACKAGE_NAME}"
            NAMESPACE "${PACKAGE_NAME}::"
            COMPONENT headers
        )

        get_property(component_names GLOBAL PROPERTY "${PACKAGE_NAME}_COMPONENTS")
        list(APPEND component_names "${COMPONENT_NAME}")
        set_property(GLOBAL PROPERTY "${PACKAGE_NAME}_COMPONENTS" ${component_names})
    endif()
endmacro()

# caribou_add_target_to_component (
#     TARGET_NAME           <target_name>
#     COMPONENT_NAME        <component_name>
#     PACKAGE_NAME          <package_name>
#     TARGET_VERSION        <version>           # defaults to ${PACKAGE_NAME}_VERSION
#     PUBLIC_HEADERS        <[file1, ...]>
#     PRIVATE_HEADERS       <[file1, ...]>
#     HEADER_SRC_PREFIX     <[relative_path]>   # Path from where the relative path to the header files will be
#                                               # computed. Default to "${CMAKE_CURRENT_SOURCE_DIR}/../..".
#     HEADER_BUILD_PREFIX   <[relative_path]>   # Path where the header file will be configured. Default to "${CMAKE_BINARY_DIR}/include/__P__"
#                                               # where __P__ is the relative path between ${HEADER_SRC_PREFIX} and the header file path.
#     HEADER_INSTALL_PREFIX <[relative_path]>   # Path where the header file will be installed. Default to "include/__P__".
#                                               # where __P__ is the relative path between ${HEADER_SRC_PREFIX} and the header file path.
#     RUNTIME_DESTINATION   <[folder_path]>     # default to "bin"
#     LIBRARY_DESTINATION   <[folder_path]>     # default to "lib"
#     ARCHIVE_DESTINATION   <[folder_path]>     # default to "lib"
# )
#
# Adds a target to a cmake component created using caribou_add_component.
# The PUBLIC_HEADERS and PRIVATE_HEADERS arguments allows to link header
# files to the build directory. If the header files extension suffix is ".in:",
# for example "config.h.in", they will be configured by cmake instead of linked.
# Only files in the PUBLIC_HEADERS will be installed.
#
# The HEADER_PREFIX can be used to specify the relative prefix path where the header files will be configured/installed.
# It is default to ${PACKAGE_NAME}/${COMPONENT_NAME}.
macro(caribou_add_target_to_component)
    set(oneValueArgs TARGET_NAME COMPONENT_NAME PACKAGE_NAME HEADER_SRC_PREFIX HEADER_BUILD_PREFIX HEADER_INSTALL_PREFIX RUNTIME_DESTINATION LIBRARY_DESTINATION ARCHIVE_DESTINATION)
    set(multiValueArgs PUBLIC_HEADERS PRIVATE_HEADERS)
    set(optionalArgs )
    cmake_parse_arguments("ARG" "${optionalArgs}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

    # Header src prefix
    if (NOT ARG_HEADER_SRC_PREFIX)
        set(ARG_HEADER_SRC_PREFIX "${CMAKE_CURRENT_SOURCE_DIR}/../..")
    endif()

    # Header build prefix
    if (NOT ARG_HEADER_BUILD_PREFIX)
        set(ARG_HEADER_BUILD_PREFIX "${CMAKE_BINARY_DIR}/include/__P__")
    endif()

    # Header install prefix
    if (NOT ARG_HEADER_INSTALL_PREFIX)
        set(ARG_HEADER_INSTALL_PREFIX "include/__P__")
    endif()

    # Target version
    set(TARGET_VERSION ${ARG_TARGET_VERSION})
    if (NOT TARGET_VERSION)
        set(TARGET_VERSION ${${ARG_PACKAGE_NAME}_VERSION})
    endif()

    # Installation paths
    if (NOT ARG_RUNTIME_DESTINATION)
        set(ARG_RUNTIME_DESTINATION "bin")
    endif()
    if (NOT ARG_LIBRARY_DESTINATION)
        set(ARG_LIBRARY_DESTINATION "lib")
    endif()
    if (NOT ARG_ARCHIVE_DESTINATION)
        set(ARG_ARCHIVE_DESTINATION "lib")
    endif()

    # The behavior is different whether or not we have an interface target (without binary) or not
    get_target_property(target_type ${ARG_TARGET_NAME} TYPE)
    set(VISIBILITY PUBLIC)
    if (target_type STREQUAL "INTERFACE_LIBRARY")
        set(TARGET_INTERFACE TRUE)
        set(VISIBILITY INTERFACE)
    endif()

    # Set the preprocessor token CARIBOU_BUILD_{TARGET_NAME}
    if (NOT TARGET_INTERFACE)
        string(TOUPPER "${ARG_TARGET_NAME}" TARGET_NAME_UPPER)
        string(REPLACE "." "_" TARGET_NAME_UPPER "${TARGET_NAME_UPPER}")
        target_compile_definitions(${ARG_TARGET_NAME} PRIVATE "-DCARIBOU_BUILD_${TARGET_NAME_UPPER}")
    endif()

    # Target properties
    if (NOT TARGET_INTERFACE)
        if (TARGET_VERSION)
            set_target_properties(${ARG_TARGET_NAME} PROPERTIES VERSION "${TARGET_VERSION}")
        endif()
        set_target_properties(${ARG_TARGET_NAME} PROPERTIES POSITION_INDEPENDENT_CODE ON)
    endif()

    # Compile features
    target_compile_features(${ARG_TARGET_NAME} ${VISIBILITY} cxx_std_17)

    # Include directories
    list(APPEND build_prefixes "${CMAKE_BINARY_DIR}/include")
    if (ARG_HEADER_SRC_PREFIX)
        list(APPEND build_prefixes "${ARG_HEADER_SRC_PREFIX}")
    endif ()
    target_include_directories(${ARG_TARGET_NAME} ${VISIBILITY} "$<BUILD_INTERFACE:${build_prefixes}>")
    target_include_directories(${ARG_TARGET_NAME} ${VISIBILITY} $<INSTALL_INTERFACE:include>)

    # Install target
    install(
        TARGETS ${ARG_TARGET_NAME}
        EXPORT  ${ARG_COMPONENT_NAME}Targets
        RUNTIME DESTINATION "${ARG_RUNTIME_DESTINATION}" COMPONENT applications
        LIBRARY DESTINATION "${ARG_LIBRARY_DESTINATION}" COMPONENT libraries
        ARCHIVE DESTINATION "${ARG_ARCHIVE_DESTINATION}" COMPONENT libraries
        PUBLIC_HEADER DESTINATION "include" COMPONENT headers
    )

    # Configure all header files (public or private)
    set(header_files ${ARG_PUBLIC_HEADERS} ${ARG_PRIVATE_HEADERS})
    foreach(header_file ${header_files})
        set(dir_from_src "")
        get_filename_component(output_filename ${header_file} NAME)

        if (ARG_HEADER_SRC_PREFIX)
            file(RELATIVE_PATH path_from_src "${ARG_HEADER_SRC_PREFIX}" "${CMAKE_CURRENT_SOURCE_DIR}/${header_file}")
            get_filename_component(dir_from_src ${path_from_src} DIRECTORY)
            get_filename_component(output_filename ${path_from_src} NAME)
        endif()

        if (output_filename MATCHES "\\.in$")
            string(REGEX REPLACE "\\.in$" "" configured_filename ${output_filename})
            string(REPLACE "__P__" "${dir_from_src}" configured_path "${ARG_HEADER_BUILD_PREFIX}")
            configure_file("${header_file}" "${configured_path}/${configured_filename}")
        endif()
    endforeach()

    # Install only public header files
    set(header_files ${ARG_PUBLIC_HEADERS})
    foreach(header_file ${header_files})
        set(dir_from_src "")
        get_filename_component(output_filename ${header_file} NAME)

        if (ARG_HEADER_SRC_PREFIX)
            file(RELATIVE_PATH path_from_src "${ARG_HEADER_SRC_PREFIX}" "${CMAKE_CURRENT_SOURCE_DIR}/${header_file}")
            get_filename_component(dir_from_src ${path_from_src} DIRECTORY)
            get_filename_component(output_filename ${path_from_src} NAME)
        endif()

        if (output_filename MATCHES "\\.in$")
            string(REGEX REPLACE "\\.in$" "" filename ${output_filename})
            string(REPLACE "__P__" "${dir_from_src}" configuration_path "${ARG_HEADER_BUILD_PREFIX}")
            string(REPLACE "__P__" "${dir_from_src}" installation_path "${ARG_HEADER_INSTALL_PREFIX}")
            install(
                FILES "${configuration_path}/${filename}"
                DESTINATION "${installation_path}"
                COMPONENT headers
            )
        else()
            string(REPLACE "__P__" "${dir_from_src}" installation_path "${ARG_HEADER_INSTALL_PREFIX}")
            install(
                FILES "${header_file}"
                DESTINATION "${installation_path}"
                COMPONENT headers
            )
        endif()
    endforeach()

endmacro()
