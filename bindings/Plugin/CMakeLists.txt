project(Bindings.Plugin)

set(HEADER_FILES
    Forcefield/FictitiousGridElasticForce.h
    Forcefield/FictitiousGridHyperelasticForcefield.h
    Forcefield/HexahedronElasticForce.h
    Forcefield/HyperelasticForcefield.h
    Ode/LegacyStaticODESolver.h
    Ode/StaticODESolver.h
    Solver/ConjugateGradientSolver.h
    Topology/FictitiousGrid.h
)

set(SOURCE_FILES
    SofaCaribou.cpp
    Forcefield/FictitiousGridElasticForce.cpp
    Forcefield/FictitiousGridHyperelasticForcefield.cpp
    Forcefield/HexahedronElasticForce.cpp
    Forcefield/HyperelasticForcefield.cpp
    Ode/LegacyStaticODESolver.cpp
    Ode/StaticODESolver.cpp
    Solver/ConjugateGradientSolver.cpp
    Topology/FictitiousGrid.cpp
)

set(PYTHON_FILES
    __init__.py
)

find_package(SofaPython3 REQUIRED)
caribou_add_python_module(SofaCaribou
                          TARGET_NAME       ${PROJECT_NAME}
                          TARGET_ALIAS      Caribou::${PROJECT_NAME}
                          TARGET_DEPENDS    Caribou::Geometry Caribou::Bindings.Config SofaCaribou SofaPython3::Bindings.Sofa.Core
                          COMPONENT_NAME    Bindings
                          PACKAGE_NAME      Caribou
                          DESTINATION       SofaCaribou
                          SOURCE_FILES      ${SOURCE_FILES}
                          PUBLIC_HEADERS    ${HEADER_FILES}
                          HEADER_PREFIX     SofaCaribou/Python
                          PYTHON_FILES      ${PYTHON_FILES}
                          PYTHON_TEST_FILES ${PYTHON_TEST_FILES}
                          )

# Avoid an annoying warning since the bindings code are in "caribou/bindings/*" but the headers are installed into
# "Caribou/Bindings/*" to be consistent with SOFA and the rest of caribou libraries.
target_compile_options(${PROJECT_NAME} PRIVATE "-Wno-nonportable-include-path")