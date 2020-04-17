set(DIRECTORY SofaCaribou)

set(PYTHON_FILES
    __init__.py
)

set(PYTHON_TEST_FILES
)

set(HEADER_FILES
    Forcefield/FictitiousGridElasticForce.h
    Forcefield/FictitiousGridHyperelasticForcefield.h
    Forcefield/HexahedronElasticForce.h
    Forcefield/HyperelasticForcefield.h
    Ode/StaticODESolver.h
    Topology/FictitiousGrid.h
)

set(SOURCE_FILES
    SofaCaribou.cpp
    Forcefield/FictitiousGridElasticForce.cpp
    Forcefield/FictitiousGridHyperelasticForcefield.cpp
    Forcefield/HexahedronElasticForce.cpp
    Forcefield/HyperelasticForcefield.cpp
    Ode/StaticODESolver.cpp
    Topology/FictitiousGrid.cpp
)

unset(PYTHON_MODULE_PREFIX)
unset(PYTHON_MODULE_EXTENSION)


if (CARIBOU_WITH_PYTHON_3)
    find_package(SofaPython3 REQUIRED)
    caribou_add_python_module(SofaCaribouPython
        PYTHON_VERSION 3
        TARGET_NAME SofaCaribouPython
        DESTINATION SofaCaribou
        SOURCE_FILES ${SOURCE_FILES}
        PYTHON_FILES ${PYTHON_FILES}
        PYTHON_TEST_FILES ${PYTHON_TEST_FILES})

    target_link_libraries(SofaCaribouPython PRIVATE SofaCaribou Caribou::Geometry SofaPython3::Plugin SofaPython3::Bindings.Sofa.Core)
    target_include_directories(SofaCaribouPython
        PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/>"
        PUBLIC $<INSTALL_INTERFACE:include>)
    set_property(TARGET SofaCaribouPython PROPERTY CXX_STANDARD 17)
endif()
