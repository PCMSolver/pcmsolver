file(READ "${CMAKE_SOURCE_DIR}/VERSION" PCMSOLVER_VERSION)
string(STRIP "${PCMSOLVER_VERSION}" PCMSOLVER_VERSION)

# reset GIT_REVISION
set(GIT_REVISION)

# if GIT_HASH exists then this is exported code
# in this case we read git hash from this file and set DEVELOPMENT_CODE to false
execute_process(COMMAND git rev-parse --abbrev-ref HEAD OUTPUT_VARIABLE is_release)
if(${is_release} STREQUAL "release")
    set(DEVELOPMENT_CODE FALSE)
else()
    set(DEVELOPMENT_CODE TRUE)
    add_definitions(-DWAVELET_DEVELOPMENT)
endif()
