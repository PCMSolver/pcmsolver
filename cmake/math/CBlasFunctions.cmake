macro(init_vendor_cblas _vendor)
if (${_vendor}_CBLAS_INCLUDE_DIRS AND 
		${_vendor}_CBLAS_LIBRARIES)
	set(${_vendor}_CBLAS_FIND_QUIETLY TRUE)
endif ()

if (NOT DEFINED ${_vendor}_ROOT)
	if (EXISTS $ENV{${_vendor}_ROOT})
		set(${_vendor}_ROOT $ENV{${_vendor}_ROOT})
	elseif (DEFINED CBLAS_ROOT)
		set(${_vendor}_ROOT ${CBLAS_ROOT})
	elseif (EXISTS $ENV{CBLAS_ROOT})
		set(${_vendor}_ROOT $ENV{CBLAS_ROOT})
	endif()
endif()
endmacro()

FUNCTION(find_cblas_include_dirs _vendor _name)
if (DEFINED ${_vendor}_ROOT)
	find_path(${_vendor}_CBLAS_INCLUDE_DIRS
		NAMES ${_name}
		PATHS ${${_vendor}_ROOT}
		PATH_SUFFIXES include
		NO_DEFAULT_PATH
		)
else()
	find_path(${_vendor}_CBLAS_INCLUDE_DIRS NAMES ${_name})
endif()	
if (${_vendor}_CBLAS_INCLUDE_DIRS)
	if (NOT DEFINED CBLAS_H_NAME)
		set(CBLAS_H_NAME ${_name} CACHE STRING "Name of CBLAS header file")
	endif()
endif()
ENDFUNCTION()

FUNCTION(find_cblas_libraries _vendor _suffix)
if (NOT ${_vendor}_CBLAS_LIBRARIES)
	foreach(_libname ${ARGN})
		if (DEFINED ${_vendor}_ROOT)
			find_library(_lib ${_libname}
				PATHS ${${_vendor}_ROOT}
				PATH_SUFFIXES ${_suffix}
				NO_DEFAULT_PATH
				)
		endif()
			find_library(_lib ${_libname} 
				PATH_SUFFIXES ${_suffix}
				)
		if(DEFINED _lib)
			set(_libs ${_libs} ${_lib})
		endif()
		unset(_lib CACHE)
	endforeach()
endif()

if (DEFINED _libs)
	set(${_vendor}_CBLAS_LIBRARIES ${_libs}
		CACHE STRING "${_vendor} BLAS libraries" FORCE)
endif()
ENDFUNCTION()
