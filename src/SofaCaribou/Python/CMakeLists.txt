project(SofaCaribou.Python)

set(HEADER_FILES
    Forcefield/HexahedronElasticForce.h
    Forcefield/HyperelasticForcefield.h
    Ode/LegacyStaticODESolver.h
    Ode/StaticODESolver.h
    Solver/ConjugateGradientSolver.h
    Topology/CaribouTopology.h
    Topology/FictitiousGrid.h
)

set(SOURCE_FILES
    SofaCaribou.cpp
    Forcefield/HexahedronElasticForce.cpp
    Forcefield/HyperelasticForcefield.cpp
    Ode/LegacyStaticODESolver.cpp
    Ode/StaticODESolver.cpp
    Solver/ConjugateGradientSolver.cpp
    Topology/CaribouTopology.cpp
    Topology/FictitiousGrid.cpp
)

set(PYTHON_FILES
    __init__.py
)

find_package(SofaPython3 REQUIRED)
caribou_add_python_module(SofaCaribou
                          TARGET_NAME       ${PROJECT_NAME}
                          TARGET_ALIAS      Caribou::Python.SofaCaribou
                          TARGET_DEPENDS    Caribou::Geometry Caribou::Python.Config SofaCaribou SofaPython3::Bindings.Sofa.Core
                          COMPONENT_NAME    SofaCaribou
                          PACKAGE_NAME      SofaCaribou
                          DESTINATION       SofaCaribou
                          SOURCE_FILES      ${SOURCE_FILES}
                          PUBLIC_HEADERS    ${HEADER_FILES}
                          PYTHON_FILES      ${PYTHON_FILES}
                          PYTHON_TEST_FILES ${PYTHON_TEST_FILES}
                          )

