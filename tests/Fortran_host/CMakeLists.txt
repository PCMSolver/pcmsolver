link_directories(${PROJECT_BINARY_DIR}/lib)
add_executable(Fortran_host Fortran_host.f90 Fortran_host-modules.f90)
if(STATIC_LIBRARY_ONLY)
  target_link_libraries(Fortran_host pcm-static)
else()
  target_link_libraries(Fortran_host pcm-shared)
endif()
include_directories(${PROJECT_SOURCE_DIR}/tests/C_host
                    ${PROJECT_SOURCE_DIR}/api
                    ${PROJECT_BINARY_DIR}/include)
set_target_properties(Fortran_host PROPERTIES LINKER_LANGUAGE Fortran)
add_test(NAME Fortran_host
         COMMAND Fortran_host
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(Fortran_host PROPERTIES LABELS "api;Fortran_host")
