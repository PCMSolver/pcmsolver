file(READ ${PROJECT_SOURCE_DIR}/README.md _readme)
string(REGEX MATCH "([0-9]+)\\.([0-9]+)\\.([0-9]+)(-[a-z,A-Z,0-9]+)?" PCMSolver_VERSION ${_readme})
string(REPLACE "." ";" _version_list ${PCMSolver_VERSION})
list(GET _version_list 0 PROJECT_VERSION_MAJOR)
list(GET _version_list 1 PROJECT_VERSION_MINOR)
list(GET _version_list 2 PROJECT_VERSION_PATCH)
# Separate on the dash
string(FIND ${PROJECT_VERSION_PATCH} "-" _has_describe)
if(NOT _has_describe EQUAL -1)
  string(REPLACE "-" ";" _repatch ${PROJECT_VERSION_PATCH})
  list(GET _repatch 1 PROJECT_VERSION_DESCRIBE)
endif()

message(STATUS "${BoldGreen}PCMSolver v${PCMSolver_VERSION}${ColourReset}")
