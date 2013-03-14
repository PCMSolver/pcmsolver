macro(config_math_service _SERVICE)
    set(ENABLE_${_SERVICE}
        ENABLE_${_SERVICE}
        CACHE BOOL
        "Enable ${_SERVICE}"
        )
    set(${_SERVICE}_FOUND FALSE)
    if(ENABLE_${_SERVICE})
        find_package(${_SERVICE})
    endif()
    if(${_SERVICE}_FOUND)
        set(LIBS
            ${LIBS}
            ${${_SERVICE}_LIBRARIES}
            )
    else()
        if(ENABLE_${_SERVICE})
            message("-- No external ${_SERVICE} library found")
        endif()
        message("-- Using own ${_SERVICE} implementation (slow)")
        add_definitions(-DUSE_OWN_${_SERVICE})
    endif()
endmacro()

set(MATH_LANG "Fortran")

set(EXPLICIT_LIBS
    "${EXPLICIT_LIBS}"
    CACHE STRING
    "User set math libraries"
    FORCE
    )

if(EXPLICIT_LIBS)
    set(LIBS
        ${LIBS}
        ${EXPLICIT_LIBS}
        )
    message("-- User set math libraries: ${EXPLICIT_LIBS}")
else()
    config_math_service(LAPACK)
    config_math_service(BLAS)
endif()
