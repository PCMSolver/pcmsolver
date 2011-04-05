# - Find a CBLAS library
#
# This module defines
#  CBLAS_INCLUDE_DIRS, where to find cblas.h (or equivalent)
#  CBLAS_LIBRARIES, the libraries to link against to use CBLAS.
#  compiling
#  CBLAS_FOUND, If false, do not try to use CBLAS.
# also defined, but not for general use are
# None of the above will be defined unless CBLAS can be found.
# 

#=============================================================================
# Copyright 2010 Jonas Juselius <jonas.juselius@uit.no>
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

if (CBLAS_INCLUDE_DIRS AND CBLAS_LIBRARIES)
  set(CBLAS_FIND_QUIETLY TRUE)
endif ()

if (DEFINED ENV{CBLAS_VENDOR}) 
	set(CBLAS_VENDOR $ENV{CBLAS_VENDOR})
endif()

if (DEFINED CBLAS_VENDOR) 
	set (_cblas_packages ${CBLAS_VENDOR})
else()
	set (_cblas_packages Mkl Goto Atlas Generic)
endif()

foreach (_cblas ${_cblas_packages})
	string(TOUPPER ${_cblas} BLASNAME)
	set(blas_package ${_cblas}CBLAS)
	string(TOUPPER ${_cblas}_CBLAS BLASLIB)

	if (EXISTS ${CBLAS_ROOT})
		if (NOT $ENV{${BLASNAME}_ROOT})
			set(ENV{${BLASNAME}_ROOT} ${CBLAS_ROOT})
		endif()
	elseif (EXISTS $ENV{CBLAS_ROOT})
		if (NOT $ENV{${BLASNAME}_ROOT})
			set(ENV{${BLASNAME}_ROOT} $ENV{CBLAS_ROOT})
		endif()
	endif()
	
	find_package(${blas_package} QUIET)
	
	if (${BLASNAME}CBLAS_FOUND)
		set(CBLAS_INCLUDE_DIRS ${${BLASLIB}_INCLUDE_DIRS}
			CACHE STRING "Current CBLAS include directories")
		set(CBLAS_LIBRARIES ${${BLASLIB}_LIBRARIES}
			CACHE STRING "Current CBLAS libraries")
		unset(${BLASLIB}_INCLUDE_DIRS CACHE)
		unset(${BLASLIB}_LIBRARIES CACHE) 
		break()
	endif()
endforeach()
unset(_cblas)
unset(_cblas_packages)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(CBLAS DEFAULT_MSG
                                  CBLAS_INCLUDE_DIRS CBLAS_LIBRARIES)
if (CBLAS_FOUND)
	set(HAVE_CBLAS ON CACHE INTERNAL "Defined if CBLAS is available")
endif()

if (NOT DEFINED CBLAS_H_NAME)
	set(CBLAS_H_NAME cblas.h CACHE STRING "Name of CBLAS header file")
endif()

mark_as_advanced(CBLAS_INCLUDE_DIRS CBLAS_LIBRARIES CBLAS_H_NAME)

unset(BLASLIB)
unset(blas_package)
unset(BLASNAME)
