include_directories(${PROJECT_SOURCE_DIR}/api)

install(FILES ${PROJECT_SOURCE_DIR}/api/pcmsolver.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PN})
install(FILES ${PROJECT_SOURCE_DIR}/api/PCMInput.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PN})

if(ENABLE_FORTRAN_API)
    add_definitions(-DENABLE_FORTRAN_API)
    add_library(fortran_bindings OBJECT ${PROJECT_SOURCE_DIR}/api/pcmsolver.F90)
    list(APPEND _objects $<TARGET_OBJECTS:fortran_bindings>)
endif()
