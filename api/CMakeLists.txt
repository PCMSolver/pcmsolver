install(FILES ${PROJECT_SOURCE_DIR}/api/pcmsolver.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME})
install(FILES ${PROJECT_SOURCE_DIR}/api/PCMInput.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME})

if(ENABLE_FORTRAN_API)
    add_library(fortran_bindings OBJECT ${PROJECT_SOURCE_DIR}/api/pcmsolver.f90)
    list(APPEND _objects $<TARGET_OBJECTS:fortran_bindings>)
endif()
