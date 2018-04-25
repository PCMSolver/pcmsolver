add_executable(C_host C_host.c)
target_include_directories(C_host
  PRIVATE
    ${PROJECT_SOURCE_DIR}/tests/C_host
  )
target_link_libraries(C_host PCMSolver)
target_link_libraries(C_host m)
set_target_properties(C_host
  PROPERTIES
    LINKER_LANGUAGE C
  )
add_test(
  NAME
    C_host
  COMMAND
    C_host
  WORKING_DIRECTORY
    ${CMAKE_CURRENT_BINARY_DIR}
  )
set_tests_properties(C_host
  PROPERTIES
    LABELS "api;C_host"
  )

add_executable(fail-C_host fail-C_host.c)
target_include_directories(fail-C_host
  PRIVATE
    ${PROJECT_SOURCE_DIR}/tests/C_host
  )
target_link_libraries(fail-C_host PCMSolver)
target_link_libraries(fail-C_host m)
set_target_properties(fail-C_host
  PROPERTIES
    LINKER_LANGUAGE C
  )
add_test(
  NAME
    fail-C_host
  COMMAND
    fail-C_host
  WORKING_DIRECTORY
    ${CMAKE_CURRENT_BINARY_DIR}
  )
set_tests_properties(fail-C_host
  PROPERTIES
    LABELS "api;C_host;fail"
    WILL_FAIL TRUE
  )
