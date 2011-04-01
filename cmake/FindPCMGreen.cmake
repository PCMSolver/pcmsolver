# - Find the PCMGreen includes and library
#
# This module defines
#  PCMGREEN_INCLUDE_DIRS, where to find GreensFunction.h, etc.
#  PCMGREEN_LIBRARIES, the libraries to link against to use XCFun.
#  PCMGREEN_DEFINITIONS - You should add_definitons(${PCMGREEN_DEFINITIONS}) before
#  compiling
#  PCMGREEN_FOUND, If false, do not try to use PCMGreen.
# also defined, but not for general use are
# None of the above will be defined unles PCMGreen can be found.
# 

#=============================================================================
# Copyright 2010 Jonas Juselius <jonas.juselius@uit.no>
# Modified by Luca Frediani from XCFUN version
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

if (PCMGREEN_INCLUDE_DIRS AND PCMGREEN_LIBRARIES)
	set(PCMGREEN_FIND_QUIETLY TRUE)
endif ()

find_path(PCMGREEN_INCLUDE_DIRS
  NAMES GreensFunction.h
  PATHS ${PCMGREEN_ROOT} $ENV{PCMGREEN_ROOT}
  PATH_SUFFIXES include
  NO_DEFAULT_PATH
)
find_path(PCMGREEN_INCLUDE_DIRS GreensFunction.h)

find_path(PCMGREEN_LIBRARIES 
  PATHS ${PCMGREEN_ROOT} $ENV{PCMGREEN_ROOT}
  PATH_SUFFIXES lib
  NO_DEFAULT_PATH
)
find_library(PCMGREEN_LIBRARIES green)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PCMGREEN DEFAULT_MSG
                                  PCMGREEN_INCLUDE_DIR PCMGREEN_LIBRARIES)

mark_as_advanced(PCMGREEN_INCLUDE_DIR PCMGREEN_LIBRARIES)
