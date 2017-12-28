install(FILES ${PROJECT_SOURCE_DIR}/api/pcmsolver.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME})
install(FILES ${PROJECT_SOURCE_DIR}/api/PCMInput.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME})
install(FILES ${PROJECT_SOURCE_DIR}/api/PCMSolverInput.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME})

if(ENABLE_Fortran_API)
    add_library(fortran_bindings OBJECT ${PROJECT_SOURCE_DIR}/api/pcmsolver.f90)
    set_target_properties(fortran_bindings
      PROPERTIES
        INCLUDE_DIRECTORIES ""
        POSITION_INDEPENDENT_CODE 1
      )
    list(APPEND _objects $<TARGET_OBJECTS:fortran_bindings>)
endif()
