add_executable(Fortran_host Fortran_host.F90 Fortran_host-modules.F90)
target_link_libraries(Fortran_host pcm getkw ${ZLIB_LIBRARIES})
set_target_properties(Fortran_host PROPERTIES INCLUDE_DIRECTORIES ${PROJECT_SOURCE_DIR}/tests/Fortran_host
                                                INCLUDE_DIRECTORIES ${PROJECT_SOURCE_DIR}/api
                                                LINKER_LANGUAGE Fortran)
add_test(NAME Fortran_host COMMAND Fortran_host)
set_tests_properties(Fortran_host PROPERTIES LABELS "api;Fortran_host")
