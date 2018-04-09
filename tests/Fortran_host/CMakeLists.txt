add_executable(Fortran_host
  Fortran_host.f90
  Fortran_host-modules.f90
  ${PROJECT_SOURCE_DIR}/api/pcmsolver.f90
  )
target_link_libraries(Fortran_host PCMSolver)
set_target_properties(Fortran_host PROPERTIES LINKER_LANGUAGE Fortran)
add_test(
  NAME
    Fortran_host
  COMMAND
    Fortran_host
  WORKING_DIRECTORY
    ${CMAKE_CURRENT_BINARY_DIR}
  )
set_tests_properties(Fortran_host
  PROPERTIES
    LABELS "api;Fortran_host"
  )
