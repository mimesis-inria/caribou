cmake_minimum_required(VERSION 3.1)
project(Caribou.Python)

set(DIRECTORY Caribou)

set(PYTHON_FILES
        __init__.py
        )

set(PYTHON_TEST_FILES
        )

set(SOURCE_FILES
        )

if (CARIBOU_WITH_PYTHON_2)
    set(PACKAGE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/site-packages/2.7")

    file(COPY ${PYTHON_FILES} DESTINATION "${PACKAGE_INSTALL_PREFIX}/${DIRECTORY}")
    foreach(t ${PYTHON_TEST_FILES})
        configure_file(${t} "${CMAKE_BINARY_DIR}/python_tests/2.7/${DIRECTORY}/${t}")
    endforeach()
    install(FILES ${PYTHON_FILES} DESTINATION site-packages/2.7/${DIRECTORY})
endif()