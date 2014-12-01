# - Find libdl
# Find the native libdl includes and library
#
#  LIBDL_INCLUDE_DIR - where to find dlfcn.h
#  LIBDL_LIBRARIES   - List of libraries when using libdl.
#  LIBDL_FOUND       - True if libdl found.

if(LIBDL_INCLUDE_DIR)
  # Already in cache, be silent
  set(LIBDL_FIND_QUIETLY TRUE)
endif(LIBDL_INCLUDE_DIR)

find_path(LIBDL_INCLUDE_DIR dlfcn.h
	  PATHS
	  /usr/include
	  /usr/local/include
	  )

set(LIBDL_NAMES dl libdl ltdl libltdl)
find_library(LIBDL_LIBRARY 
	     NAMES ${LIBDL_NAMES}
	     PATHS
	     /usr/lib
	     /usr/local/lib
	     /lib
	     )

# handle the QUIETLY and REQUIRED arguments and set LIBDL_FOUND to TRUE if 
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(libdl DEFAULT_MSG LIBDL_LIBRARY LIBDL_INCLUDE_DIR)

if(LIBDL_FOUND)
  set(LIBDL_LIBRARIES ${LIBDL_LIBRARY})
else(LIBDL_FOUND)
  set(LIBDL_LIBRARIES)
endif(LIBDL_FOUND)

mark_as_advanced(LIBDL_LIBRARY LIBDL_INCLUDE_DIR)
