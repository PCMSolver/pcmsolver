function(get_Python_version_string pyVersion) 
    execute_process(COMMAND "${PYTHON_EXECUTABLE}" -c
                            "import sys; sys.stdout.write(';'.join([str(x) for x in sys.version_info[:3]]))"
                    OUTPUT_VARIABLE _VERSION
                    RESULT_VARIABLE _PYTHON_VERSION_RESULT
                    ERROR_QUIET)
    if(NOT _PYTHON_VERSION_RESULT)
        string(REPLACE ";" "." pyVersion "${_VERSION}")
        list(GET _VERSION 0 PYTHON_VERSION_MAJOR)
        list(GET _VERSION 1 PYTHON_VERSION_MINOR)
        list(GET _VERSION 2 PYTHON_VERSION_PATCH)
        if(PYTHON_VERSION_PATCH EQUAL 0)
            # it's called "Python 2.7", not "2.7.0"
            string(REGEX REPLACE "\\.0$" "" pyVersion "${pyVersion}")
        endif()
    else()
        # sys.version predates sys.version_info, so use that
        execute_process(COMMAND "${PYTHON_EXECUTABLE}" -c "import sys; sys.stdout.write(sys.version)"
                        OUTPUT_VARIABLE _VERSION
                        RESULT_VARIABLE _PYTHON_VERSION_RESULT
                        ERROR_QUIET)
        if(NOT _PYTHON_VERSION_RESULT)
            string(REGEX REPLACE " .*" "" pyVersion "${_VERSION}")
            string(REGEX REPLACE "^([0-9]+)\\.[0-9]+.*" "\\1" PYTHON_VERSION_MAJOR "${pyVersion}")
            string(REGEX REPLACE "^[0-9]+\\.([0-9])+.*" "\\1" PYTHON_VERSION_MINOR "${pyVersion}")
            if(pyVersion MATCHES "^[0-9]+\\.[0-9]+\\.[0-9]+.*")
                string(REGEX REPLACE "^[0-9]+\\.[0-9]+\\.([0-9]+).*" "\\1" PYTHON_VERSION_PATCH "${pyVersion}")
            else()
                set(PYTHON_VERSION_PATCH "0")
            endif()
        else()
            # sys.version was first documented for Python 1.5, so assume
            # this is older.
            set(pyVersion "1.4")
            set(PYTHON_VERSION_MAJOR "1")
            set(PYTHON_VERSION_MAJOR "4")
            set(PYTHON_VERSION_MAJOR "0")
        endif()
    endif()
    set(pyVersion "${pyVersion}" PARENT_SCOPE)
    unset(_PYTHON_VERSION_RESULT)
    unset(_VERSION)
endfunction()

function(check_Python_compiles pyCompiles)
   set(_bindir  "${PROJECT_BINARY_DIR}/pyCompiles")	
   set(_srcfile "${PROJECT_SOURCE_DIR}/cmake/embedded_Python.cpp")
   
   try_compile(pyCompiles "${_bindir}" "${_srcfile}" 
             CMAKE_FLAGS 
	     "-DINCLUDE_DIRECTORIES=${PYTHON_INCLUDE_DIRS}"
             "-DLINK_LIBRARIES=${PYTHON_LIBRARIES}")

   set(pyCompiles "${pyCompiles}" PARENT_SCOPE)	
endfunction()
