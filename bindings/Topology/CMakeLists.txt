project(Bindings.Topology)

set(PYTHON_FILES
    __init__.py
)

set(PYTHON_TEST_FILES
    pytest/Caribou.Topology.Test.py
)

set(HEADER_FILES
    BarycentricContainer.h
    Domain.h
)

set(SOURCE_FILES
    Grid.cpp
    HashGrid.cpp
    Mesh.cpp
    IO/IO.cpp
    Topology.cpp
)


caribou_add_python_module(Topology
                          TARGET_NAME       ${PROJECT_NAME}
                          TARGET_ALIAS      Caribou::${PROJECT_NAME}
                          TARGET_DEPENDS    Caribou::Topology Caribou::Bindings.Config
                          COMPONENT_NAME    Bindings
                          PACKAGE_NAME      Caribou
                          DESTINATION       Caribou/Topology
                          SOURCE_FILES      ${SOURCE_FILES}
                          PUBLIC_HEADERS    ${HEADER_FILES}
                          HEADER_PREFIX     Caribou/Bindings/Topology
                          PYTHON_FILES      ${PYTHON_FILES}
                          PYTHON_TEST_FILES ${PYTHON_TEST_FILES}
                          )