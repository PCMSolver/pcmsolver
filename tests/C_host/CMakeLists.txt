link_directories(${PROJECT_BINARY_DIR}/lib)
add_executable(C_host C_host.c)
if(STATIC_LIBRARY_ONLY)
  target_link_libraries(C_host pcm-static ${ZLIB_LIBRARIES})
else()
  target_link_libraries(C_host pcm-shared)
endif()
target_link_libraries(C_host m)
set_target_properties(C_host PROPERTIES INCLUDE_DIRECTORIES ${PROJECT_SOURCE_DIR}/tests/C_host
                                        INCLUDE_DIRECTORIES ${PROJECT_SOURCE_DIR}/api
                                        LINKER_LANGUAGE C)
add_test(NAME C_host COMMAND C_host)
set_tests_properties(C_host PROPERTIES LABELS "api;C_host")
