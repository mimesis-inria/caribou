project(Caribou.Python.Topology)

set(PYTHON_FILES
    __init__.py
)

set(PYTHON_TEST_FILES
    pytest/Caribou_Topology_Test.py
)

set(HEADER_FILES
    BarycentricContainer.h
    Domain.h
)

set(SOURCE_FILES
    Domain.cpp
    Grid.cpp
    HashGrid.cpp
    Mesh.cpp
    IO/IO.cpp
    Topology.cpp
)


caribou_add_python_module(Topology
                          TARGET_NAME       ${PROJECT_NAME}
                          TARGET_ALIAS      Caribou::Python.Topology
                          TARGET_DEPENDS    Caribou::Topology Caribou::Python.Config
                          COMPONENT_NAME    Python
                          PACKAGE_NAME      Caribou
                          DESTINATION       Caribou/Topology
                          SOURCE_FILES      ${SOURCE_FILES}
                          PUBLIC_HEADERS    ${HEADER_FILES}
                          HEADER_SRC_PREFIX "${CMAKE_CURRENT_SOURCE_DIR}/../../.."
                          PYTHON_FILES      ${PYTHON_FILES}
                          PYTHON_TEST_FILES ${PYTHON_TEST_FILES}
                          )

