add_executable(C_host.x C_host.c)
target_link_libraries(C_host.x pcm getkw ${ZLIB_LIBRARIES})
set_target_properties(C_host.x PROPERTIES INCLUDE_DIRECTORIES ${PROJECT_SOURCE_DIR}/tests/C_host
                                          INCLUDE_DIRECTORIES ${PROJECT_SOURCE_DIR}/api
                                          LINKER_LANGUAGE C)
add_test(NAME C_host COMMAND C_host.x)
set_tests_properties(C_host PROPERTIES LABELS "api;C_host")
