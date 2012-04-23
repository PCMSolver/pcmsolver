set(MATH_LANG "Fortran")

set(BLAS_FOUND   FALSE)
set(LAPACK_FOUND FALSE)

if(NOT ENABLE_INTERNAL_MATH)

    # user sets
    set(USERDEFINED_MATH
        "${USERDEFINED_MATH}"
        CACHE STRING
        "User set math libraries"
        FORCE
        )
    if(NOT "${USERDEFINED_MATH}" STREQUAL "")
        set(MATH_LIBS
            "${USERDEFINED_MATH}"
            CACHE STRING
            "User set math libraries"
            FORCE
            )
        message("-- User set math libraries: ${MATH_LIBS}")
        set(BLAS_FOUND   TRUE)
        set(LAPACK_FOUND TRUE)
        set(LIBS
            ${LIBS}
            ${MATH_LIBS}
            )
    endif()

    # try to find the best libraries using environment variables
    if(NOT BLAS_FOUND)
        find_package(BLAS)
        if(BLAS_FOUND)
            set(LIBS
                ${LIBS}
                ${BLAS_LIBRARIES}
                )
        endif()
    endif()
    if(NOT LAPACK_FOUND)
        find_package(LAPACK)
        if(LAPACK_FOUND)
            set(LIBS
                ${LIBS}
                ${LAPACK_LIBRARIES}
                )
        endif()
    endif()

    if(NOT BLAS_FOUND)
        message("-- No external BLAS library found")
    endif()
    if(NOT LAPACK_FOUND)
        message("-- No external LAPACK library found")
    endif()

endif()
