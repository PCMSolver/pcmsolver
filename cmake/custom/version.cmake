add_custom_target(update_version
  ALL
  COMMAND
    ${PYTHON_EXECUTABLE} tools/versioner.py --metaout ${CMAKE_CURRENT_BINARY_DIR}/metadata.py
                                            --cmakeout ${CMAKE_CURRENT_BINARY_DIR}/metadata.cmake
  COMMAND
    ${CMAKE_COMMAND} -DWTO="${CMAKE_CURRENT_BINARY_DIR}/${CMAKECONFIG_INSTALL_DIR}"
          -DPN="PCMSolver"
          -P ${CMAKE_CURRENT_BINARY_DIR}/metadata.cmake
  WORKING_DIRECTORY
    ${CMAKE_CURRENT_SOURCE_DIR}
  COMMENT
    "Generating version info"
  )
#install(FILES ${CMAKE_CURRENT_BINARY_DIR}/metadata.py
#        DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}${PYMOD_INSTALL_LIBDIR}/pylibefp)

file(READ ${PROJECT_SOURCE_DIR}/README.md _readme)
string(REGEX MATCH "([0-9]+)\\.([0-9]+)\\.([0-9]+)(-[a-z,A-Z,0-9]+)?" PCMSolver_VERSION ${_readme})
string(REPLACE "." ";" _version_list ${PCMSolver_VERSION})
list(GET _version_list 0 PROJECT_VERSION_MAJOR)
list(GET _version_list 1 PROJECT_VERSION_MINOR)
list(GET _version_list 2 _patch_describe)
# Get PROJECT_VERSION_PATCH
string(FIND ${_patch_describe} "-" _has_describe)
if(_has_describe GREATER -1)
  string(REGEX REPLACE "-[a-z,A-Z,0-9]+" "" PROJECT_VERSION_PATCH ${_patch_describe})
else()
  set(PROJECT_VERSION_PATCH ${_patch_describe})
endif()

message(STATUS "${BoldGreen}PCMSolver v${PCMSolver_VERSION}${ColourReset}")
