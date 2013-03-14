# - Find a LAPACK library
#
# This module will first look in LAPACK_ROOT before considering the default
# system pahts.
# The linker language can be defined by setting the varable LAPACK_LANG
#
# This module defines:
#
#  LAPACK_INCLUDE_DIRS Where to find lapack.h (or equivalent)
#  LAPACK_LIBRARIES Libraries to link against to use LAPACK
#  LAPACK_FOUND Defined if LAPACK is available
#  HAVE_LAPACK To be used in #ifdefs
#
# None of the above will be defined unless LAPACK can be found.
#
#=============================================================================
# Copyright 2011 Jonas Juselius <jonas.juselius@uit.no>
#                Radovan Bast   <radovan.bast@uit.no>
#
# Distributed under the OSI-approved BSD License (the "License");
# see accompanying file Copyright.txt for details.
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.
#=============================================================================
# (To distributed this file outside of CMake, substitute the full
#  License text for the above reference.)

include(MathLibFunctions)

if (EXISTS $ENV{MATH_ROOT})
    if (NOT DEFINED LAPACK_ROOT})
        set(LAPACK_ROOT $ENV{MATH_ROOT})
    endif()
endif()

if (EXISTS $ENV{LAPACK_ROOT})
    if (NOT DEFINED LAPACK_ROOT})
        set(LAPACK_ROOT $ENV{LAPACK_ROOT})
    endif()
endif()

# BLAS and LAPACK often go together
if (NOT DEFINED LAPACK_ROOT})
    if (DEFINED BLAS_ROOT})
        set(LAPACK_ROOT ${BLAS_ROOT})
    elseif (EXISTS $ENV{BLAS_ROOT})
        set(LAPACK_ROOT $ENV{BLAS_ROOT})
    endif()
endif()

# Default names for the headers

if (LAPACK_INCLUDE_DIRS AND LAPACK_LIBRARIES)
  set(LAPACK_FIND_QUIETLY TRUE)
endif ()

if (NOT LAPACK_FIND_COMPONENTS)
    if (DEFINED LAPACK_TYPE)
        set(LAPACK_FIND_COMPONENTS ${LAPACK_TYPE})
    elseif(ENABLE_64BIT_INTEGERS)
        set(LAPACK_FIND_COMPONENTS MKL)
    else()
        set(LAPACK_FIND_COMPONENTS MKL Atlas ACML default)
    endif()
endif()

function(find_lapack)
    foreach (lapack ${LAPACK_FIND_COMPONENTS})
        if (${lapack} STREQUAL "MKL")
            find_mkl()
        elseif (${lapack} STREQUAL "Atlas")
            find_atlas()
        elseif (${lapack} STREQUAL "ACML")
            find_acml()
        else()
            find_generic()
        endif()
        if (LAPACK_FOUND)
            break()
        endif()
    endforeach()
endfunction()

macro(find_generic)
    if (MATH_LANG STREQUAL "C")
        find_math_header(lapack clapack.h)
    endif()
    find_math_libs(lapack lapack)
    cache_math_result(lapack generic)
endmacro()

macro(find_acml)
    set(MATH_LIBRARY_PATH_SUFFIXES libso)
    set(MATH_INCLUDE_PATH_SUFFIXES)
    if (MATH_LANG STREQUAL "C")
        find_math_header(lapack clapack.h)
    endif()
    find_math_libs(lapack acml)
    cache_math_result(lapack ACML)
endmacro()

macro(find_atlas)
    set(MATH_LIBRARY_PATH_SUFFIXES 
        atlas atlas-base atlas-base/atlas)
    set(MATH_INCLUDE_PATH_SUFFIXES atlas)

    if (MATH_LANG STREQUAL "C")
        find_math_header(lapack clapack.h)
        find_math_libs(lapack lapack_atlas lapack)
    endif()
    find_math_libs(lapack lapack_atlas lapack)
    cache_math_result(lapack Atlas) 
endmacro()

macro(find_mkl)
    if (MATH_LANG STREQUAL "C")
        find_math_header(lapack mkl_clapack.h)
    endif()

    if(${CMAKE_HOST_SYSTEM_PROCESSOR} STREQUAL "x86_64")
        set(MATH_LIBRARY_PATH_SUFFIXES intel64 em64t)
    else()
        set(MATH_LIBRARY_PATH_SUFFIXES ia32 32)
    endif()

    if(${CMAKE_HOST_SYSTEM_PROCESSOR} STREQUAL "x86_64")
        if(ENABLE_64BIT_INTEGERS)
            set(lapack_libs mkl_lapack95_ilp64)
        else()
            set(lapack_libs mkl_lapack95_lp64)
        endif()
    else()
        set(lapack_libs mkl_lapack95)
    endif()
    find_math_libs(lapack ${lapack_libs})

    if (NOT LAPACK_LIBRARIES)
        set(lapack_libs mkl_lapack)
        find_math_libs(lapack ${lapack_libs})
    endif()

    if (LAPACK_LIBRARIES)
        set(LAPACK_LIBRARIES 
            -Wl,--start-group ${LAPACK_LIBRARIES} -Wl,--end-group)
    endif()
    cache_math_result(lapack MKL)
    unset(lapack_libs)
endmacro()

find_lapack()
if(LAPACK_FOUND)
   find_package_message(LAPACK "Found LAPACK: ${LAPACK_TYPE}"
       "[${LAPACK_LIBRARIES}]")
endif()

