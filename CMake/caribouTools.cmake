function(caribou_add_python_module NAME)
    set(options QUIET)
    set(oneValueArgs PREFIX TESTS_PREFIX DESTINATION TESTS_DESTINATION PYTHON_VERSION TARGET_NAME)
    set(multiValueArgs SOURCE_FILES PYTHON_FILES PYTHON_TEST_FILES DEPENDS)

    cmake_parse_arguments(A "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

    if (NOT DEFINED A_PYTHON_VERSION)
        message(FATAL_ERROR "A python version must be provided.")
    endif()

    if (A_TARGET_NAME)
        set(TARGET_NAME ${A_TARGET_NAME})
    else()
        set(TARGET_NAME "${NAME}.python.${A_PYTHON_VERSION}")
    endif()


    if (A_TESTS_PREFIX)
        set(TESTS_PREFIX "${A_TESTS_PREFIX}")
    else()
        set(TESTS_PREFIX "bin/python_tests/python${A_PYTHON_VERSION}")
    endif()

    if (A_DESTINATION)
        set(DESTINATION "${A_DESTINATION}")
    else()
        set(DESTINATION "")
    endif ()

    if (A_TESTS_DESTINATION)
        set(TESTS_DESTINATION "${A_TESTS_DESTINATION}")
    else()
        set(TESTS_DESTINATION "")
    endif ()

    if (A_SOURCE_FILES)
        set(PYBIND11_PYTHON_VERSION ${A_PYTHON_VERSION})
        set(PYBIND11_CPP_STANDARD -std=c++17)
        unset(PYTHON_MODULE_PREFIX)
        unset(PYTHON_MODULE_EXTENSION)
        unset(PYTHON_EXECUTABLE)
        unset(PYTHONLIBS_FOUND)
        unset(PYTHON_MODULE_EXTENSION)
        unset(PythonLibsNew_FIND_REQUIRED )

        project(${TARGET_NAME})

        if (${A_PYTHON_VERSION} VERSION_LESS "3.0")
            find_package(PythonLibs ${A_PYTHON_VERSION} QUIET)
        else()
            find_package(Python3 QUIET COMPONENTS Development)
            set(PYTHON_EXECUTABLE ${Python3_EXECUTABLE})
        endif ()

        if (A_PREFIX)
            set(PREFIX "${A_PREFIX}")
        else()
            if ("${Python3_VERSION_MINOR}" STREQUAL "0")
                set(PREFIX "python${Python3_VERSION_MAJOR}/site-packages")
            else()
                set(PREFIX "python${Python3_VERSION_MAJOR}.${Python3_VERSION_MINOR}/site-packages")
            endif()
        endif()

        find_package(pybind11 CONFIG REQUIRED)

        set(MODULENAME ${NAME})

        if (NOT A_QUIET)
            message(STATUS "${NAME} with python ${A_PYTHON_VERSION} support (python version ${PYTHON_VERSION_STRING}, pybind11 version ${pybind11_VERSION})")
        endif ()

        pybind11_add_module(${TARGET_NAME} SHARED "${A_SOURCE_FILES}")

        target_link_libraries(${TARGET_NAME} PUBLIC ${A_DEPENDS} ${PYTHON_LIBRARIES} pybind11::module)

        target_include_directories(${TARGET_NAME}
                                   PUBLIC $<INSTALL_INTERFACE:include>
                                   )
        if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
            target_compile_options(${TARGET_NAME} PRIVATE -fsized-deallocation)
        endif()

        target_compile_options(${TARGET_NAME} PRIVATE -Dregister=)
        set_property(TARGET ${TARGET_NAME} PROPERTY CXX_STANDARD 17)

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

    if (A_PYTHON_FILES)
        foreach(t ${A_PYTHON_FILES})
            configure_file(${t} "${CMAKE_BINARY_DIR}/lib/${PREFIX}/${DESTINATION}/${t}")
            install(FILES "${CMAKE_BINARY_DIR}/lib/${PREFIX}/${DESTINATION}/${t}" DESTINATION lib/${PREFIX}/${DESTINATION})
        endforeach()
    endif()

    if (A_PYTHON_TEST_FILES)
        foreach(t ${A_PYTHON_TEST_FILES})
            configure_file(${t} "${CMAKE_BINARY_DIR}/${TESTS_PREFIX}/${TESTS_DESTINATION}/${t}")
            install(FILES "${CMAKE_BINARY_DIR}/${TESTS_PREFIX}/${TESTS_DESTINATION}/${t}" DESTINATION ${TESTS_PREFIX}/${TESTS_DESTINATION})
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
