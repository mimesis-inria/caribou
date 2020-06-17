set(DIRECTORY Caribou/Topology)

set(PYTHON_FILES
        __init__.py
)

set(PYTHON_TEST_FILES
        test/Caribou.Topology.Test.py
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

unset(PYTHON_MODULE_PREFIX)
unset(PYTHON_MODULE_EXTENSION)

if (CARIBOU_WITH_PYTHON_3)
    caribou_add_python_module(CaribouTopologyPython
                              DESTINATION ${DIRECTORY}
                              SOURCE_FILES ${SOURCE_FILES}
                              PYTHON_FILES ${PYTHON_FILES}
                              PYTHON_TEST_FILES ${PYTHON_TEST_FILES}
                              DEPENDS Topology Geometry)
endif()
