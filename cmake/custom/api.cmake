include_directories(${PROJECT_SOURCE_DIR}/api)

file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/include)

file(COPY ${PROJECT_SOURCE_DIR}/api/pcmsolver.h DESTINATION ${PROJECT_BINARY_DIR}/include)

if(ENABLE_FORTRAN_API)
    set_property(GLOBAL APPEND PROPERTY PCMSolver_Fortran_SOURCES ${PROJECT_SOURCE_DIR}/api/pcmsolver.F90)
endif()
