message("-- Fortran compiler is  \"${CMAKE_Fortran_COMPILER}\" (\"${CMAKE_Fortran_COMPILER_ID}\") ")
message("--       C compiler is  \"${CMAKE_C_COMPILER}\" (\"${CMAKE_C_COMPILER_ID}\") ")

# Fortran compilers

if(CMAKE_Fortran_COMPILER_ID MATCHES GNU) # this is gfortran
    set(CMAKE_Fortran_FLAGS         "-DVAR_GFORTRAN -DGFORTRAN=445")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -fbacktrace")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -ffast-math -funroll-loops -ftree-vectorize")
    if(ENABLE_64BIT_INTEGERS)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -fdefault-integer-8"
            )
    endif()
    if(ENABLE_BOUNDS_CHECK)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -fbounds-check"
            )
    endif()
    if(ENABLE_CODE_COVERAGE)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -ftest-coverage"
            )
    endif()
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES G95)
    set(CMAKE_Fortran_FLAGS         "-Wno=155 -fno-second-underscore -DVAR_G95 -fsloppy-char")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -ftrace=full")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
    if(ENABLE_64BIT_INTEGERS)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -i8"
            )
    endif()
    if(ENABLE_BOUNDS_CHECK)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -Wall -fbounds-check"
            )
    endif()
    if(ENABLE_CODE_COVERAGE)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS}"
            )
    endif()
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
    add_definitions(-DVAR_IFORT)
    set(CMAKE_Fortran_FLAGS         "-g -w -fpp -assume byterecl")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-O0")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -xW -ip")
    if(ENABLE_64BIT_INTEGERS)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -i8"
            )
    endif()
    if(ENABLE_BOUNDS_CHECK)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -check bounds -fpstkchk -check pointers -check uninit -check output_conversion -traceback"
            )
    endif()
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES PGI)
    set(CMAKE_Fortran_FLAGS         "-DVAR_PGF90")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-g -O0 -Mframe")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -mcmodel=medium -fast -Munroll")
    if(ENABLE_64BIT_INTEGERS)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -m64 -i8"
            )
    else()
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -m32"
            )
    endif()
    if(ENABLE_BOUNDS_CHECK)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} "
            )
    endif()
    if(ENABLE_CODE_COVERAGE)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} "
            )
    endif()
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES XL)
    set(CMAKE_Fortran_FLAGS         "-qzerosize -qextname")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-g")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
    if(ENABLE_64BIT_INTEGERS)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -q64"
            )
    endif()

    set_source_files_properties(${FREE_FORTRAN_SOURCES}
        PROPERTIES COMPILE_FLAGS
        "-qfree"
        )
    set_source_files_properties(${FIXED_FORTRAN_SOURCES}
        PROPERTIES COMPILE_FLAGS
        "-qfixed"
        )
endif()

# C compilers

if(CMAKE_C_COMPILER_ID MATCHES GNU)
    set(CMAKE_C_FLAGS         "-std=c99 -DRESTRICT=restrict -DFUNDERSCORE=1")
    set(CMAKE_C_FLAGS_DEBUG   "-O0 -g")
    set(CMAKE_C_FLAGS_RELEASE "-O3 -ffast-math -funroll-loops -ftree-vectorize -Wno-unused")
endif()

if(CMAKE_C_COMPILER_ID MATCHES Intel)
    set(CMAKE_C_FLAGS         "-g -wd981 -wd279 -wd383 -vec-report0 -wd1572 -wd177 -restrict -DRESTRICT=restrict")
    set(CMAKE_C_FLAGS_DEBUG   "-O0")
    set(CMAKE_C_FLAGS_RELEASE "-O3 -ip")
    set(CMAKE_C_LINK_FLAGS "${CMAKE_C_LINK_FLAGS} -shared-intel")
endif()

if(CMAKE_C_COMPILER_ID MATCHES PGI)
    set(CMAKE_C_FLAGS         " ")
    set(CMAKE_C_FLAGS_DEBUG   "-g -O0")
    set(CMAKE_C_FLAGS_RELEASE "-O3 -fast -Munroll -Mvect=idiom -c9x -DRESTRICT=restrict")
endif()

if(CMAKE_C_COMPILER_ID MATCHES XL)
    set(CMAKE_C_FLAGS         " ")
    set(CMAKE_C_FLAGS_DEBUG   " ")
    set(CMAKE_C_FLAGS_RELEASE " ")
endif()
