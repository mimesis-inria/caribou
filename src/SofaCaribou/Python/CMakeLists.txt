set(DIRECTORY SofaCaribou)

set(PYTHON_FILES
    __init__.py
    )

set(PYTHON_TEST_FILES
    )

set(HEADER_FILES
    Forcefield/FictitiousGridElasticForce.h
    Forcefield/HexahedronElasticForce.h
    Topology/FictitiousGrid.h
    )

set(SOURCE_FILES
    SofaCaribou.cpp
    Forcefield/FictitiousGridElasticForce.cpp
    Forcefield/HexahedronElasticForce.cpp
    Topology/FictitiousGrid.cpp
    )

unset(PYTHON_MODULE_PREFIX)
unset(PYTHON_MODULE_EXTENSION)


if (CARIBOU_WITH_PYTHON_3)
    caribou_add_python_module(SofaCaribouPython
        PYTHON_VERSION 3
        TARGET_NAME SofaCaribouPython
        DESTINATION SofaCaribou
        SOURCE_FILES ${SOURCE_FILES}
        PYTHON_FILES ${PYTHON_FILES}
        PYTHON_TEST_FILES ${PYTHON_TEST_FILES})

    target_link_libraries(SofaCaribouPython PRIVATE SofaCaribou Caribou::Geometry)
    target_include_directories(SofaCaribouPython
        PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/>"
        PUBLIC $<INSTALL_INTERFACE:include>)
    set_property(TARGET SofaCaribouPython PROPERTY CXX_STANDARD 17)
endif()
