# - Find a CBLAS library
#
# This module defines
#  CBLAS_INCLUDE_DIRS, where to find cblas.h (or equivalent)
#  CBLAS_LIBRARIES, the libraries to link against to use CBLAS.
#  compiling
#  CBLAS_FOUND, If false, do not try to use CBLAS.
# also defined, but not for general use are
# None of the above will be defined unles GENERIC_CBLAS can be found.
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

include(CBlasFunctions)
init_vendor_cblas(GENERIC)
find_cblas_include_dirs(GENERIC cblas.h)
find_cblas_libraries(GENERIC lib cblas)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GENERIC_CBLAS DEFAULT_MSG
	GENERIC_CBLAS_INCLUDE_DIRS GENERIC_CBLAS_LIBRARIES)

if (GENERIC_CBLAS_FOUND)
	set(CBLAS_FOUND TRUE)
	set(CBLAS_INCLUDE_DIRS ${GENERIC_CBLAS_INCLUDE_DIRS})
	set(CBLAS_LIBRARIES ${GENERIC_CBLAS_LIBRARIES})
endif()

mark_as_advanced(GENERIC_CBLAS_INCLUDE_DIRS GENERIC_CBLAS_LIBRARIES)
