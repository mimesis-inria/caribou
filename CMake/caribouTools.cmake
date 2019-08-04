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

    if (A_PREFIX)
        set(PREFIX "${A_PREFIX}")
    else()
        set(PREFIX "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/site-packages/${A_PYTHON_VERSION}")
    endif()

    if (A_TESTS_PREFIX)
        set(TESTS_PREFIX "${A_TESTS_PREFIX}")
    else()
        set(TESTS_PREFIX "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/python_tests/${A_PYTHON_VERSION}")
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

        find_package(pybind11 CONFIG REQUIRED)

        set(MODULENAME ${NAME})

        if (NOT A_QUIET)
            message(STATUS "${NAME} with python ${A_PYTHON_VERSION} support (python version ${PYTHON_VERSION_STRING}, pybind11 version ${pybind11_VERSION})")
        endif ()

        pybind11_add_module(${TARGET_NAME} SHARED "${A_SOURCE_FILES}")

        target_link_libraries(${TARGET_NAME} PUBLIC ${A_DEPENDS} ${PYTHON_LIBRARIES} pybind11::module)

        target_include_directories(${TARGET_NAME}
                                   PUBLIC "$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src/>"
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
            LIBRARY_OUTPUT_DIRECTORY "${PREFIX}/${DESTINATION}"
        )
    endif ()

    if (A_PYTHON_FILES)
        foreach(t ${A_PYTHON_FILES})
            configure_file(${t} "${PREFIX}/${DESTINATION}/${t}")
            install(FILES "${PREFIX}/${DESTINATION}/${t}" DESTINATION ${PREFIX}/${DESTINATION})
        endforeach()
    endif()

    if (A_PYTHON_TEST_FILES)
        foreach(t ${A_PYTHON_TEST_FILES})
            configure_file(${t} "${TESTS_PREFIX}/${TESTS_DESTINATION}/${t}")
            install(FILES "${TESTS_PREFIX}/${TESTS_DESTINATION}/${t}" DESTINATION ${TESTS_PREFIX}/${TESTS_DESTINATION})
        endforeach()
    endif()

endfunction()