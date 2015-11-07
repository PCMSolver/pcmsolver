add_executable(Fortran_host.x Fortran_host.F90 Fortran_host-modules.F90)
target_link_libraries(Fortran_host.x pcm getkw ${ZLIB_LIBRARIES})
set_target_properties(Fortran_host.x PROPERTIES INCLUDE_DIRECTORIES ${PROJECT_SOURCE_DIR}/tests/Fortran_host
                                                INCLUDE_DIRECTORIES ${PROJECT_SOURCE_DIR}/api
                                                LINKER_LANGUAGE Fortran)
add_test(NAME Fortran_host COMMAND Fortran_host.x)
set_tests_properties(Fortran_host PROPERTIES LABELS "api;Fortran_host")
