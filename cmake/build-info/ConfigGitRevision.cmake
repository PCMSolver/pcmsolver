# only in the development code we get the hash from git
# in the released code we read it from file, in this case
# it is already set at this stage
if(DEVELOPMENT_CODE)
    find_package(Git)
    if(GIT_FOUND)
        execute_process(
            COMMAND ${GIT_EXECUTABLE} rev-list --max-count=1 --abbrev-commit HEAD
            OUTPUT_VARIABLE GIT_REVISION
            ERROR_QUIET
        )

        if(GIT_REVISION)
            string(STRIP ${GIT_REVISION} GIT_REVISION)
        endif()

        execute_process(
            COMMAND ${GIT_EXECUTABLE} log --max-count=1 HEAD
            OUTPUT_VARIABLE GIT_LAST_COMMIT
            ERROR_QUIET
        )

        if(GIT_LAST_COMMIT)
            execute_process(
                COMMAND ${GIT_EXECUTABLE} --no-pager log -1 --pretty=format:%an HEAD
                OUTPUT_VARIABLE GIT_LAST_COMMIT_AUTHOR
                ERROR_QUIET
                )
            execute_process(
                COMMAND ${GIT_EXECUTABLE} --no-pager log -1 --pretty=format:%ae HEAD
                OUTPUT_VARIABLE GIT_LAST_COMMIT_AUTHOR_EMAIL
                ERROR_QUIET
                )
            execute_process(
                COMMAND ${GIT_EXECUTABLE} --no-pager log -1 --pretty=format:%ad HEAD
                OUTPUT_VARIABLE GIT_LAST_COMMIT_DATE
                ERROR_QUIET
                )
        endif()

        execute_process(
            COMMAND ${PYTHON_EXECUTABLE} -c "import subprocess; import re; print(re.search(r'\\*.*', subprocess.Popen(['${GIT_EXECUTABLE}', 'branch'], stdout=subprocess.PIPE).communicate()[0], re.MULTILINE).group())"
            OUTPUT_VARIABLE GIT_BRANCH
            ERROR_QUIET
        )

        if(GIT_BRANCH)
            string(REPLACE "*" "" GIT_BRANCH ${GIT_BRANCH})
            string(STRIP ${GIT_BRANCH} GIT_BRANCH)
        endif()

        execute_process(
            COMMAND ${GIT_EXECUTABLE} status -uno
            OUTPUT_FILE ${PROJECT_BINARY_DIR}/GIT_STATUS_AT_BUILD
            ERROR_QUIET
        )
    endif()
endif()
