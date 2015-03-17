
if(DEVELOPMENT_CODE)
    include(FindGit)
    if(GIT_FOUND)
        add_custom_target(
            git_update
            COMMAND ${GIT_EXECUTABLE} submodule init
            COMMAND ${GIT_EXECUTABLE} submodule update
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
    else()
            message("-- Git not found. You need Git for the Git submodule mechanism to work.")
    endif()
endif()

include(ExternalProject)

macro(add_external _project)

    if(DEVELOPMENT_CODE AND GIT_FOUND)
        set(UPDATE_COMMAND ${GIT_EXECUTABLE} submodule update)
    else()
        set(UPDATE_COMMAND echo)
    endif()

    add_custom_target(
        check_external_timestamp_${_project}
	COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_SOURCE_DIR}/cmake/check_external_timestamp.py
                       ${CMAKE_BINARY_DIR}/external/${_project}-stamp/${_project}-configure
                       ${CMAKE_BINARY_DIR}/external/${_project}-stamp
                       ${CMAKE_SOURCE_DIR}/external/${_project}
    )

    # Set the testing command
    set(_testing_command "${ARGV1}")
    # Set the external project args.
    # In case no argument was passed for this, use the global variable
    set(_external_project_cmake_args "")
    if("${ARGV2}" STREQUAL "")
	    set(_external_project_cmake_args "${ExternalProjectCMakeArgs}")
    else()
	    set(_external_project_cmake_args "${ARGV2}")
    endif()

    if("${_testing_command}" STREQUAL "")
       ExternalProject_Add(${_project}
           DOWNLOAD_COMMAND ${UPDATE_COMMAND}
           DOWNLOAD_DIR ${CMAKE_SOURCE_DIR}
           SOURCE_DIR ${CMAKE_SOURCE_DIR}/external/${_project}
           BINARY_DIR ${CMAKE_BINARY_DIR}/external/${_project}-build
           STAMP_DIR ${CMAKE_BINARY_DIR}/external/${_project}-stamp
           TMP_DIR ${CMAKE_BINARY_DIR}/external/${_project}-tmp
           INSTALL_DIR ${CMAKE_BINARY_DIR}/external
           CMAKE_ARGS ${_external_project_cmake_args}
           )
    else()
       # For unfathomable reasons, CMake expects the TEST_COMMAND to be a ;-separated list...
       separate_arguments(_testing_command)
       ExternalProject_Add(${_project}
           DOWNLOAD_COMMAND ${UPDATE_COMMAND}
           DOWNLOAD_DIR ${CMAKE_SOURCE_DIR}
           SOURCE_DIR ${CMAKE_SOURCE_DIR}/external/${_project}
           BINARY_DIR ${CMAKE_BINARY_DIR}/external/${_project}-build
           STAMP_DIR ${CMAKE_BINARY_DIR}/external/${_project}-stamp
           TMP_DIR ${CMAKE_BINARY_DIR}/external/${_project}-tmp
           INSTALL_DIR ${CMAKE_BINARY_DIR}/external
           CMAKE_ARGS ${_external_project_cmake_args}
	   TEST_BEFORE_INSTALL 1
	   TEST_COMMAND "${_testing_command}"
	   LOG_TEST 1
           )
    endif()

    include_directories(${CMAKE_BINARY_DIR}/external/${_project}-build)
    link_directories(${CMAKE_BINARY_DIR}/external/lib)

    if(DEVELOPMENT_CODE)
        add_dependencies(${_project} git_update)
        add_dependencies(check_external_timestamp_${_project} git_update)
    endif()
    add_dependencies(${_project} check_external_timestamp_${_project})
endmacro()
