set(DIRECTORY Caribou/Topology)

set(PYTHON_FILES
    __init__.py
)

set(PYTHON_TEST_FILES
)

set(HEADER_FILES
        Grid.h
        HashGrid.h
)

set(SOURCE_FILES
        Grid.cpp
        HashGrid.cpp
        Topology.cpp
)

unset(PYTHON_MODULE_PREFIX)
unset(PYTHON_MODULE_EXTENSION)

if (CARIBOU_WITH_PYTHON_2)
    caribou_add_python_module(CaribouTopologyPython
                              PYTHON_VERSION 2.7
                              DESTINATION ${DIRECTORY}
                              SOURCE_FILES ${SOURCE_FILES}
                              PYTHON_FILES ${PYTHON_FILES}
                              PYTHON_TEST_FILES ${PYTHON_TEST_FILES}
                              DEPENDS Topology Geometry)
endif()

if (CARIBOU_WITH_PYTHON_3)
    caribou_add_python_module(CaribouTopologyPython
                              PYTHON_VERSION 3
                              DESTINATION ${DIRECTORY}
                              SOURCE_FILES ${SOURCE_FILES}
                              PYTHON_FILES ${PYTHON_FILES}
                              PYTHON_TEST_FILES ${PYTHON_TEST_FILES}
                              DEPENDS Topology Geometry)
endif()
