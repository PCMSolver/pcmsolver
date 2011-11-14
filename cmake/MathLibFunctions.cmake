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

if (NOT MATH_LANG)
	set(MATH_LANG C)
elseif(MATH_LANG STREQUAL "C" OR MATH_LANG STREQUAL "CXX")
	set(MATH_LANG C)
elseif(NOT MATH_LANG STREQUAL "Fortran")
	message(FATAL_ERROR "Invalid math library linker language: ${MATH_LANG}")
endif()

macro(find_math_header _service)
	string(TOUPPER ${_service} _SERVICE)
	if (${_service}_h)
		find_path(${_service}_include_dirs
			NAMES ${${_service}_h}
			PATHS ${${_SERVICE}_ROOT}
			PATH_SUFFIXES include ${path_suffixes}
			NO_DEFAULT_PATH
			)
		find_path(${_service}_include_dirs
			NAMES ${${_service}_h}
			PATH_SUFFIXES include
			)
	endif()
endmacro()

macro(find_math_libs _service)
	string(TOUPPER ${_service} _SERVICE)
	foreach(${_service}lib ${${_service}_libs})
		find_library(_lib ${${_service}lib}
			PATHS ${${_SERVICE}_ROOT}
			PATH_SUFFIXES ${path_suffixes}
			NO_DEFAULT_PATH
			)
		find_library(_lib ${${_service}lib}
			PATH_SUFFIXES ${path_suffixes}
			)
		if(_lib)
			set(${_service}_libraries ${${_service}_libraries} ${_lib})
			unset(_lib CACHE)
		else()
			break()
		endif()
	endforeach()
	unset(${_service}lib)
	unset(_lib CACHE)
endmacro()

macro(cache_math_result math_type _service)
	string(TOUPPER ${_service} _SERVICE)
	if (${_service}_h)
		set(${_SERVICE}_H ${${_service}_h})
	endif()
	if (${_service}_include_dirs)
		set(${_SERVICE}_INCLUDE_DIRS ${${_service}_include_dirs})
	endif()
	if (${_service}_libraries)
		set(${_SERVICE}_LIBRARIES ${${_service}_libraries})
	endif()
	unset(${_service}_h)
	unset(${_service}_include_dirs)
	unset(${_service}_libraries)

    set(${_SERVICE}_FIND_QUIETLY TRUE)
	if (${_SERVICE}_H)
		find_package_handle_standard_args(${_SERVICE}
			"Could NOT find ${math_type} ${_SERVICE}"
			${_SERVICE}_INCLUDE_DIRS ${_SERVICE}_LIBRARIES ${_SERVICE}_H)
	else()
		find_package_handle_standard_args(${_SERVICE}
			"Could NOT find ${math_type} ${_SERVICE}" ${_SERVICE}_LIBRARIES)
	endif()

	if (${_SERVICE}_FOUND)
		message("-- ${math_type} ${_SERVICE} found")
		set(HAVE_${_SERVICE} ON CACHE INTERNAL "Defined if ${_SERVICE} is available")
		set(${_SERVICE}_LIBRARIES ${${_SERVICE}_LIBRARIES} CACHE STRING "${_SERVICE} libraries")
		mark_as_advanced(${_SERVICE}_LIBRARIES)
		if (${_SERVICE}_H)
			set(${_SERVICE}_H ${${_SERVICE}_H} CACHE STRING "Name of ${_SERVICE} header")
			set(${_SERVICE}_INCLUDE_DIRS ${${_SERVICE}_INCLUDE_DIRS}
				CACHE STRING "${_SERVICE} include directory")
			mark_as_advanced(${_SERVICE}_INCLUDE_DIRS ${_SERVICE}_H)
		endif()
	endif()
endmacro()
