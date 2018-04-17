cmake_minimum_required(VERSION 2.8.12)
project(caribou)

set(HEADER_FILES
        src/Event/IterativeSolverEvent.h
        src/Solver/NewtonRaphsonSolver.h
)

set(SOURCE_FILES
        src/Event/IterativeSolverEvent.cpp
        src/Solver/NewtonRaphsonSolver.cpp
)

file(GLOB SCENE_FILES
        "scenes/*.py"
        "scenes/*.scn"
        "scenes/*/*.py"
        "scenes/*/*.scn"
)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic")

find_package(SofaGeneral REQUIRED)
find_package(SofaBase REQUIRED)

OPTION(WITH_PYTHON_BINDING "Compile python bindings" ON)
if (WITH_PYTHON_BINDING)
    find_package(SofaPython REQUIRED)

    set(HEADER_FILES ${HEADER_FILES}
            src/Python/PythonEventBinder.h
            src/Python/SofaCaribouBindings.h
            src/Python/Bindings/PythonEventBinder.h
    )

    set(SOURCE_FILES ${SOURCE_FILES}
            src/InitCaribou.cpp
            src/Python/PythonEventBinder.cpp
            src/Python/SofaCaribouBindings.cpp
            src/Python/Bindings/PythonEventBinder.cpp
    )

    add_definitions(-DWITH_PYTHON_BINDING)
endif ()


add_library(${PROJECT_NAME} SHARED ${HEADER_FILES} ${SOURCE_FILES} ${SCENE_FILES})
target_link_libraries(${PROJECT_NAME} SofaComponentGeneral SofaBaseTopology)

if (WITH_PYTHON_BINDING)
    target_link_libraries(${PROJECT_NAME} SofaPython)
endif ()
