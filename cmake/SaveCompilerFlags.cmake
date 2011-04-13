# Take care of updating the cache for fresh configurations

macro(SaveCompilerFlags _lang)

if (NOT DEFINED HAVE_${_lang}_FLAGS)
	mark_as_advanced(HAVE_${_lang}_FLAGS)
	
    set(CMAKE_${_lang}_FLAGS "${CMAKE_${_lang}_FLAGS}"
        CACHE STRING 
		"Flags used by the compiler during all builds." FORCE)

    set(CMAKE_${_lang}_FLAGS_DEBUG "${CMAKE_${_lang}_FLAGS_DEBUG}"
        CACHE STRING 
		"Flags used by the compiler during debug builds." FORCE)

    set(CMAKE_${_lang}_FLAGS_RELEASE "${CMAKE_${_lang}_FLAGS_RELEASE}"
        CACHE STRING 
		"Flags used by the compiler during release builds." FORCE)

	set (HAVE_${_lang}_FLAGS ON 
		CACHE BOOL
		"Flag that the default ${_lang} compiler flags have been set.")
endif()
endmacro()

