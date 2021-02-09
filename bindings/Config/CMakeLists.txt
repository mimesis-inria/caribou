project(Bindings.Config)

set(PYTHON_FILES
    __init__.py
)

set(PYTHON_TEST_FILES
)

set(HEADER_FILES
    Caribou.h
)

set(SOURCE_FILES
    Caribou.cpp
)

caribou_add_python_module(Caribou
                          TARGET_NAME       ${PROJECT_NAME}
                          TARGET_ALIAS      Caribou::${PROJECT_NAME}
                          TARGET_DEPENDS    Caribou::Config
                          COMPONENT_NAME    Bindings
                          PACKAGE_NAME      Caribou
                          DESTINATION       Caribou
                          SOURCE_FILES      ${SOURCE_FILES}
                          PUBLIC_HEADERS    ${HEADER_FILES}
                          HEADER_PREFIX     Caribou/Bindings
                          PYTHON_FILES      ${PYTHON_FILES}
)

# Avoid an annoying warning since the bindings code are in "caribou/bindings/*" but the headers are installed into
# "Caribou/Bindings/*" to be consistent with SOFA and the rest of caribou libraries.
target_compile_options(${PROJECT_NAME} PRIVATE "-Wno-nonportable-include-path")