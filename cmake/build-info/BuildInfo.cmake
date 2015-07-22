include(ConfigGitRevision)

# Generate the BuildInfo.hpp header only if the logger was enabled
if(ENABLE_LOGGER)
   execute_process(
       COMMAND ${PYTHON_EXECUTABLE} -c "import getpass; print getpass.getuser()"
       TIMEOUT 1
       OUTPUT_VARIABLE USER_NAME
       OUTPUT_STRIP_TRAILING_WHITESPACE
       )

   execute_process(
       COMMAND ${PYTHON_EXECUTABLE} -c "from socket import gethostname; print gethostname()"
       TIMEOUT 1
       OUTPUT_VARIABLE HOST_NAME
       OUTPUT_STRIP_TRAILING_WHITESPACE
       )

   configure_file(
       ${CMAKE_SOURCE_DIR}/cmake/build-info/build_info.py.in
       ${CMAKE_BINARY_DIR}/build_info.py
       )

   add_custom_target(
       generate_build_info
       COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_BINARY_DIR}/build_info.py > ${CMAKE_BINARY_DIR}/include/BuildInfo.hpp
       COMMAND ${CMAKE_COMMAND} -E remove -f ${CMAKE_BINARY_DIR}/build_info.py
       )

   add_dependencies(pcm generate_build_info)

   set_source_files_properties(${CMAKE_SOURCE_DIR}/src/utils/BuildInfo.hpp PROPERTIES GENERATED 1)
endif()
