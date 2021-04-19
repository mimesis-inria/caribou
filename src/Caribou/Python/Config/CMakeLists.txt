project(Caribou.Python.Config)

set(PYTHON_FILES
    __init__.py
)

set(PYTHON_TEST_FILES
)

set(HEADER_FILES
    Caribou.h.in
)

set(SOURCE_FILES
    Caribou.cpp
)

caribou_add_python_module(Caribou
                          TARGET_NAME       ${PROJECT_NAME}
                          TARGET_ALIAS      Caribou::Python.Config
                          TARGET_DEPENDS    Caribou::Config
                          COMPONENT_NAME    Python
                          PACKAGE_NAME      Caribou
                          DESTINATION       Caribou
                          SOURCE_FILES      ${SOURCE_FILES}
                          PUBLIC_HEADERS    ${HEADER_FILES}
                          HEADER_SRC_PREFIX "${CMAKE_CURRENT_SOURCE_DIR}/../../.."
                          HEADER_BUILD_PREFIX   "${CMAKE_BINARY_DIR}/include/Caribou/Python"
                          HEADER_INSTALL_PREFIX "include/Caribou/Python"
                          PYTHON_FILES      ${PYTHON_FILES}
)

if (NOT WIN32)
    # Avoid an annoying warning since the bindings code are in "caribou/Python/*" but the headers are installed into
    # "Caribou/Python/*" to be consistent with SOFA and the rest of caribou libraries.
    target_compile_options(${PROJECT_NAME} PRIVATE "-Wno-nonportable-include-path")
endif()