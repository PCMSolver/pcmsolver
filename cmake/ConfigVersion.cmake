set(PROJECT_VERSION_MAJOR 0)
set(PROJECT_VERSION_MINOR 2)
set(PROJECT_VERSION_PATCH 0)

# reset GIT_REVISION
set(GIT_REVISION)

# if RELEASE exists then this is exported code
# in this case set DEVELOPMENT_CODE to false
if(EXISTS ${PROJECT_SOURCE_DIR}/RELEASE)
    set(DEVELOPMENT_CODE FALSE)
else()
    set(DEVELOPMENT_CODE TRUE)
    add_definitions(-DWAVELET_DEVELOPMENT)
    add_definitions(-DTSLESS_DEVELOPMENT)
endif()
