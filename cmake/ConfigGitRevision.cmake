set (GIT_REVISION)
find_package(Git)
if (GIT_FOUND)
    execute_process(
        COMMAND ${GIT_EXECUTABLE} rev-list --abbrev-commit --max-count=1 HEAD
        OUTPUT_VARIABLE GIT_REVISION
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
        )
    string(STRIP
        ${GIT_REVISION}
        GIT_REVISION)
endif()

