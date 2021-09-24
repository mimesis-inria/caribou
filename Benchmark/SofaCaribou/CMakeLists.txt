project(SofaCaribou.Benchmark)

set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)

set(HEADER_FILES
    Forcefield/RecomputeOrStore/HyperelasticForcefieldRecomputeF.h
    Forcefield/RecomputeOrStore/HyperelasticForcefieldStoreF.h
    Forcefield/RecomputeOrStore/HyperelasticForcefieldStoreFAndS.h
)

set(TEMPLATE_FILES
    Forcefield/RecomputeOrStore/HyperelasticForcefieldRecomputeF.inl
    Forcefield/RecomputeOrStore/HyperelasticForcefieldStoreF.inl
    Forcefield/RecomputeOrStore/HyperelasticForcefieldStoreFAndS.inl
)

set(SOURCE_FILES
    Forcefield/RecomputeOrStore/RecomputeOrStore.cpp
    init.cpp
)

add_library(${PROJECT_NAME} SHARED ${SOURCE_FILES} ${TEMPLATE_FILES} ${HEADER_FILES})

target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_17)
set_target_properties(${PROJECT_NAME} PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_compile_options(${PROJECT_NAME} PRIVATE -march=native)
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

target_link_libraries(${PROJECT_NAME} SofaCaribou)
