include(FindGit)

if (GIT_FOUND)
    add_custom_target(
        git_update
        COMMAND ${GIT_EXECUTABLE} submodule init
        COMMAND ${GIT_EXECUTABLE} submodule update
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    )
else()
    if(DEVELOPMENT_CODE)
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

    ExternalProject_Add(${_project}
        DOWNLOAD_COMMAND ${UPDATE_COMMAND}
        DOWNLOAD_DIR ${PROJECT_SOURCE_DIR}
        PREFIX ${PROJECT_SOURCE_DIR}/external
        SOURCE_DIR ${PROJECT_SOURCE_DIR}/external/${_project}
        BINARY_DIR ${PROJECT_BINARY_DIR}/external/${_project}-build
        STAMP_DIR ${PROJECT_BINARY_DIR}/external/${_project}-stamp
        TMP_DIR ${PROJECT_BINARY_DIR}/external/${_project}-tmp
        INSTALL_DIR ${PROJECT_BINARY_DIR}/external
        CMAKE_ARGS ${ExternalProjectCMakeArgs}
        )
    include_directories(${PROJECT_BINARY_DIR}/external/${_project}-build)
    include_directories(${PROJECT_BINARY_DIR}/external/${_project}-build/modules)
    link_directories(${PROJECT_BINARY_DIR}/external/lib)
    link_directories(${PROJECT_BINARY_DIR}/external/${_project}-build/external/lib)
    if(DEVELOPMENT_CODE)
        add_dependencies(${_project} git_update)
    endif()

  # # remove stamps for external builds so that they are rebuilt every time
  # add_custom_command(
  #     TARGET ${_project}
  #     PRE_BUILD
  #     COMMAND rm -rf ${PROJECT_BINARY_DIR}/external/${_project}-stamp
  #     WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
  #     )
endmacro()
