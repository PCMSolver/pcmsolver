include_directories(${PROJECT_SOURCE_DIR}/api)

install(FILES ${PROJECT_SOURCE_DIR}/api/pcmsolver.h DESTINATION include)

if(ENABLE_FORTRAN_API)
    add_definitions(-DENABLE_FORTRAN_API)
    set_property(GLOBAL APPEND PROPERTY PCMSolver_Fortran_SOURCES ${PROJECT_SOURCE_DIR}/api/pcmsolver.F90)
else()
    install(FILES ${PROJECT_SOURCE_DIR}/api/PCMInput.h DESTINATION include)
endif()
