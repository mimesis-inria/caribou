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

function(caribou_add_python_module NAME)
    set(options QUIET)
    set(oneValueArgs PREFIX DESTINATION TARGET_NAME)
    set(multiValueArgs SOURCE_FILES HEADERSHEADERS PYTHON_FILES PYTHON_TEST_FILES DEPENDS)

    cmake_parse_arguments(A "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

    set(PYBIND11_PYTHON_VERSION 3)
    find_package(pybind11 QUIET CONFIG REQUIRED)

    if (A_TARGET_NAME)
        set(TARGET_NAME ${A_TARGET_NAME})
    else()
        set(TARGET_NAME "${NAME}.python.${PYTHON_LIBRARY_SUFFIX}")
    endif()

    if (A_DESTINATION)
        set(DESTINATION "${A_DESTINATION}")
    else()
        set(DESTINATION "")
    endif ()

    if (A_PREFIX)
        set(PREFIX "${A_PREFIX}")
    else()
        set(PREFIX "python${PYTHON_LIBRARY_SUFFIX}/site-packages")
    endif()

    # Fetch the current path relative to /*/src
    string(REGEX MATCH "(.*)/src" path_to_src "${CMAKE_CURRENT_SOURCE_DIR}")

    if (A_SOURCE_FILES)
        set(PYBIND11_CPP_STANDARD -std=c++17)

        project(${TARGET_NAME})

        set(MODULENAME ${NAME})

        pybind11_add_module(${TARGET_NAME} SHARED "${A_SOURCE_FILES}")
        target_link_libraries(${TARGET_NAME} PUBLIC ${A_DEPENDS} pybind11::module)
        target_include_directories(${TARGET_NAME} PUBLIC "$<BUILD_INTERFACE:${path_to_src}/>")
        target_include_directories(${TARGET_NAME} PUBLIC $<INSTALL_INTERFACE:include>)

        if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
            target_compile_options(${TARGET_NAME} PRIVATE -fsized-deallocation)
        endif()

        target_compile_options(${TARGET_NAME} PRIVATE -Dregister=)
        target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_17)

        set_target_properties(
            ${TARGET_NAME}
            PROPERTIES
            OUTPUT_NAME ${NAME}
            PREFIX "${PYTHON_MODULE_PREFIX}"
            SUFFIX "${PYTHON_MODULE_EXTENSION}"
            LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib/${PREFIX}/${DESTINATION}"
        )
        install(
               TARGETS ${TARGET_NAME}
               LIBRARY DESTINATION "lib/${PREFIX}/${DESTINATION}" COMPONENT libraries
        )
    endif ()

    foreach(header ${A_HEADERS})
        file(RELATIVE_PATH path_from_package "${path_to_src}" "${header}")
        get_filename_component(dir_from_package ${path_from_package} DIRECTORY)

        install(
            FILES
                "${header}"
            DESTINATION
                "include/${dir_from_package}"
            COMPONENT headers
        )
    endforeach()

    if (A_PYTHON_FILES)
        foreach(t ${A_PYTHON_FILES})
            configure_file(${t} "${CMAKE_BINARY_DIR}/lib/${PREFIX}/${DESTINATION}/${t}")
            install(FILES "${CMAKE_BINARY_DIR}/lib/${PREFIX}/${DESTINATION}/${t}" DESTINATION lib/${PREFIX}/${DESTINATION})
        endforeach()
    endif()

    if (A_PYTHON_TEST_FILES)
        set(CARIBOU_PYTHON_LIB_PATH "${CMAKE_BINARY_DIR}/lib/${PREFIX}/${DESTINATION}")
        string(REGEX MATCH "(.*)/site-packages" CARIBOU_PYTHON_LIB_PATH "${CARIBOU_PYTHON_LIB_PATH}")
        foreach(t ${A_PYTHON_TEST_FILES})
            get_filename_component(test_filename "${CMAKE_CURRENT_SOURCE_DIR}/${t}" NAME)
            configure_file(${t} "${CMAKE_BINARY_DIR}/bin/${test_filename}")
            install(FILES "${CMAKE_BINARY_DIR}/bin/${test_filename}" DESTINATION bin)
        endforeach()
    endif()

endfunction()


# caribou_install_target(target[, headers])
macro(caribou_install_target package target)
    set(target_version ${${target}_VERSION})
    set(package_version ${${package}_VERSION})

    if(NOT target_version VERSION_GREATER "0.0")
        set(target_version ${package_version})
    endif()


    install(TARGETS ${target}
        EXPORT ${package}Targets
        RUNTIME DESTINATION "bin" COMPONENT applications
        LIBRARY DESTINATION "lib" COMPONENT libraries
        ARCHIVE DESTINATION "lib" COMPONENT libraries
        PUBLIC_HEADER DESTINATION "include/${package}/${target}/${include_install_dir}" COMPONENT headers
    )

    get_target_property(target_type ${target} TYPE)
    if (NOT target_type STREQUAL "INTERFACE_LIBRARY")
        if(target_version VERSION_GREATER "0.0")
            set_target_properties(${target} PROPERTIES VERSION "${target_version}")
        endif()
    endif()

    if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${target}Config.cmake.in")
        configure_package_config_file(
                "${CMAKE_CURRENT_SOURCE_DIR}/${target}Config.cmake.in"
                "${CMAKE_BINARY_DIR}/cmake/${package}/${target}Config.cmake"
                INSTALL_DESTINATION
                lib/cmake/${PROJECT_NAME}
        )

        write_basic_package_version_file(
                "${CMAKE_BINARY_DIR}/cmake/${package}/${target}ConfigVersion.cmake"
                VERSION ${target_version}
                COMPATIBILITY ExactVersion)

        install(
                FILES
                "${CMAKE_BINARY_DIR}/cmake/${package}/${target}Config.cmake"
                "${CMAKE_BINARY_DIR}/cmake/${package}/${target}ConfigVersion.cmake"
                DESTINATION
                "lib/cmake/${package}"
        )
    endif()

    set(optional_headers ${ARGN})
    if (optional_headers)
        foreach(header ${optional_headers})
            file(RELATIVE_PATH path_from_package "${CMAKE_CURRENT_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}/${header}")
            get_filename_component(dir_from_package ${path_from_package} DIRECTORY)
            get_filename_component(output_filename ${path_from_package} NAME)

            if (output_filename MATCHES "\\.in$")
                string(REGEX REPLACE "\\.in$" "" output_filename ${output_filename})
                configure_file("${header}" "${CMAKE_CURRENT_BINARY_DIR}/${dir_from_package}/${output_filename}")

                install(
                    FILES "${CMAKE_CURRENT_BINARY_DIR}/${dir_from_package}/${output_filename}"
                    DESTINATION "include/${package}/${target}/${dir_from_package}"
                    COMPONENT headers
                )
            else()
                install(
                        FILES "${CMAKE_CURRENT_SOURCE_DIR}/${header}"
                        DESTINATION "include/${package}/${target}/${dir_from_package}"
                        COMPONENT headers
                )
            endif()
        endforeach()
    endif()
endmacro()
