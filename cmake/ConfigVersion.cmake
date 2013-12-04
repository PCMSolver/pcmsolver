file(READ "${CMAKE_SOURCE_DIR}/VERSION" PCMSOLVER_VERSION)
string(STRIP "${PCMSOLVER_VERSION}" PCMSOLVER_VERSION)

# reset GIT_REVISION
set(GIT_REVISION)

# if GIT_HASH exists then this is exported code
# in this case we read git hash from this file and set DEVELOPMENT_CODE to false
if(EXISTS ${CMAKE_SOURCE_DIRECTORY}/RELEASE)
    set(DEVELOPMENT_CODE FALSE)
else()
    set(DEVELOPMENT_CODE TRUE)
    add_definitions(-DWAVELET_DEVELOPMENT)
endif()
