# - Find libutil
# Find the native libutil includes and library
#
#  LIBUTIL_INCLUDE_DIR - where to find pty.h and utmp.h (UNIX) or util.h (Mac OS X)
#  LIBUTIL_LIBRARIES   - List of libraries when using libutil.
#  LIBUTIL_FOUND       - True if libutil found.

if(LIBUTIL_INCLUDE_DIR)
  # Already in cache, be silent
  set(LIBUTIL_FIND_QUIETLY TRUE)
endif(LIBUTIL_INCLUDE_DIR)

if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
   find_path(LIBUTIL_INCLUDE_DIR util.h
   	  PATHS
   	  /usr/include
   	  /usr/local/include
   	  )
else()
# Should probably be generalized for Windows
   find_path(LIBUTIL_INCLUDE_DIR pty.h
   	  PATHS
   	  /usr/include
   	  /usr/local/include
   	  )
endif()

set(LIBUTIL_NAMES util libutil)
find_library(LIBUTIL_LIBRARY 
	     NAMES ${LIBUTIL_NAMES}
	     PATHS 
	     /usr/lib
	     /usr/local/lib
	     /lib
	     )

# handle the QUIETLY and REQUIRED arguments and set LIBUTIL_FOUND to TRUE if 
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(libutil DEFAULT_MSG LIBUTIL_LIBRARY LIBUTIL_INCLUDE_DIR)

if(LIBUTIL_FOUND)
   set(LIBUTIL_LIBRARIES ${LIBUTIL_LIBRARY})
else(LIBUTIL_FOUND)
   set(LIBUTIL_LIBRARIES)
endif(LIBUTIL_FOUND)

mark_as_advanced(LIBUTIL_LIBRARY LIBUTIL_INCLUDE_DIR)
