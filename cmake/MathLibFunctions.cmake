#=============================================================================
# Copyright 2011 Jonas Juselius <jonas.juselius@uit.no>
#
# Conributions by Radovan Bast <radovan.bast@uit.no>
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

include(FindPackageHandleStandardArgs)
include(FindPackageMessage)

if (NOT MATH_LANG)
    set(MATH_LANG C)
elseif(MATH_LANG STREQUAL "C" OR MATH_LANG STREQUAL "CXX")
    set(MATH_LANG C)
elseif(NOT MATH_LANG STREQUAL "Fortran")
    message(FATAL_ERROR "Invalid math library linker language: ${MATH_LANG}")
endif()

macro(find_math_header _service _header)
    string(TOUPPER ${_service} _SERVICE)
    find_path(${_SERVICE}_INCLUDE_DIRS
        NAMES ${_header}
        PATHS ${${_SERVICE}_ROOT}
        HINTS ${${_SERVICE}_ROOT}/include
        PATH_SUFFIXES ${MATH_INCLUDE_PATH_SUFFIXES}
        NO_DEFAULT_PATH
        )
    find_path(${_SERVICE}_INCLUDE_DIRS
        NAMES ${_header}
        PATH_SUFFIXES include
        )
    set (${_SERVICE}_H ${_header})
    unset (_SERVICE)
endmacro()

macro(find_math_libs _service)
    string(TOUPPER ${_service} _SERVICE)
    if(${_SERVICE}_FOUND)
        return()
    endif()
    set(_lib)
    set(_libs)
    foreach(l ${ARGN})
        find_library(_lib
            NAMES ${l}
            PATHS ${${_SERVICE}_ROOT}
            HINTS ${${_SERVICE}_ROOT}/lib64 ${${_SERVICE}_ROOT}/lib
            PATH_SUFFIXES ${MATH_LIBRARY_PATH_SUFFIXES}
            NO_DEFAULT_PATH
            )
        find_library(_lib
            NAMES ${l}
            PATH_SUFFIXES ${MATH_LIBRARY_PATH_SUFFIXES}
            )
        if(_lib)
            set(_libs ${_libs} ${_lib})
        else()
            set(${_SERVICE}_LIBRARIES ${_SERVICE}_LIBRARIES-NOTFOUND)
            set(_libs ${_SERVICE}_LIBRARIES-NOTFOUND)
            break()
        endif()
        unset (_lib CACHE)
    endforeach()
    set(${_SERVICE}_LIBRARIES ${_libs})
    unset (_lib CACHE)
    unset (_libs CACHE)
    unset (_SERVICE)
    unset (l)
endmacro()

macro(cache_math_result _service MATH_TYPE)
    string(TOUPPER ${_service} _SERVICE)
    set(${_SERVICE}_FIND_QUIETLY TRUE)
    if (DEFINED ${_SERVICE}_INCLUDE_DIRS)
        find_package_handle_standard_args(${_SERVICE}
            "Could NOT find ${MATH_TYPE} ${_SERVICE}"
            ${_SERVICE}_LIBRARIES ${_SERVICE}_INCLUDE_DIRS)
    else()
        find_package_handle_standard_args(${_SERVICE}
            "Could NOT find ${MATH_TYPE} ${_SERVICE}" ${_SERVICE}_LIBRARIES)
    endif()

    if (${_SERVICE}_FOUND)
        set(${_SERVICE}_TYPE ${MATH_TYPE} CACHE STRING
            "${_SERVICE} type")
        mark_as_advanced(${_SERVICE}_TYPE)

        add_definitions(-DHAVE_${MATH_TYPE}_${_SERVICE})
        set(HAVE_${_SERVICE} ON CACHE INTERNAL
            "Defined if ${_SERVICE} is available"
            )
        set(${_SERVICE}_LIBRARIES ${${_SERVICE}_LIBRARIES} CACHE STRING
            "${_SERVICE} libraries"
            )
        mark_as_advanced(${_SERVICE}_LIBRARIES)
        if (DEFINED ${_SERVICE}_INCLUDE_DIRS)
            set(${_SERVICE}_H ${${_SERVICE}_H} CACHE STRING
                "${_SERVICE} header file")
            mark_as_advanced(${_SERVICE}_H)
            set(${_SERVICE}_INCLUDE_DIRS ${${_SERVICE}_INCLUDE_DIRS}
                CACHE STRING "${_SERVICE} include directory"
                )
            mark_as_advanced(${_SERVICE}_INCLUDE_DIRS)
        endif()
    else()
        set (${_SERVICE}_LIBRARIES ${_SERVICE}_LIBRARIES-NOTFOUND)
        if (DEFINED ${_SERVICE}_H)
            set (${_SERVICE}_INCLUDE_DIRS ${_SERVICE}_INCLUDE_DIRS-NOTFOUND)
            unset (${_SERVICE}_H)
        endif()
    endif()
    set(${_SERVICE}_FOUND ${${_SERVICE}_FOUND} PARENT_SCOPE)
    unset (MATH_TYPE)
    unset (_SERVICE)
endmacro()
