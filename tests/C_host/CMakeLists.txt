link_directories(${PROJECT_BINARY_DIR}/lib)
add_executable(C_host C_host.c)
add_executable(fail-C_host fail-C_host.c)
if(STATIC_LIBRARY_ONLY)
  target_link_libraries(C_host pcm-static)
  target_link_libraries(fail-C_host pcm-static)
else()
  target_link_libraries(C_host pcm-shared)
  target_link_libraries(fail-C_host pcm-shared)
endif()

include_directories(${PROJECT_SOURCE_DIR}/tests/C_host
                    ${PROJECT_SOURCE_DIR}/api
                    ${PROJECT_BINARY_DIR}/include)
target_link_libraries(C_host m)
set_target_properties(C_host PROPERTIES LINKER_LANGUAGE C)
add_test(NAME C_host COMMAND C_host
                     WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(C_host PROPERTIES LABELS "api;C_host")

target_link_libraries(fail-C_host m)
set_target_properties(fail-C_host PROPERTIES LINKER_LANGUAGE C)
add_test(NAME fail-C_host COMMAND fail-C_host
                     WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(fail-C_host PROPERTIES LABELS "api;C_host;fail")
set_tests_properties(fail-C_host PROPERTIES WILL_FAIL TRUE)
