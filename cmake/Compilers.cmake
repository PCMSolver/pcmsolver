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

if (NOT DEFINED HAVE_CXX_FLAGS)
if (CMAKE_COMPILER_IS_GNUCXX)
	set (CMAKE_CXX_FLAGS "-Wall -Wno-unknown-pragmas -Wno-sign-compare -Woverloaded-virtual -Wwrite-strings -Wno-unused")
	set (CMAKE_CXX_FLAGS_DEBUG "-O0 -g3 -DDEBUG")
	set (CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -Wno-unused")
	if (ENABLE_CODE_COVERAGE)
		set (CMAKE_CXX_FLAGS 
			"${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
		set (CMAKE_CXX_LINK_FLAGS "-fprofile-arcs -ftest-coverage")
	endif()
elseif (CMAKE_CXX_COMPILER_ID MATCHES Intel)
	set (CMAKE_CXX_FLAGS "-Wno-unknown-pragmas")
	set (CMAKE_CXX_FLAGS_DEBUG "-O0 -debug -DDEBUG")
	set (CMAKE_CXX_FLAGS_RELEASE "-debug -O3 -DNDEBUG")
	set (CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS} -shared-intel")
endif ()
SaveCompilerFlags(CXX)
endif ()

if (NOT DEFINED HAVE_C_FLAGS)
if (CMAKE_COMPILER_IS_GNUCC)
	set (CMAKE_C_FLAGS "-Wall -Wno-sign-compare")
	set (CMAKE_C_FLAGS_DEBUG "-O0 -g3 -DDEBUG")
	set (CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG -Wno-unused")
	if (ENABLE_CODE_COVERAGE)
		set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} --coverage")
	endif()
elseif (CMAKE_C_COMPILER_ID MATCHES Intel)
	set (CMAKE_C_FLAGS "")
	set (CMAKE_C_FLAGS_DEBUG "-O0 -g -DDEBUG")
	set (CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG")
	set (CMAKE_C_LINK_FLAGS "${CMAKE_C_LINK_FLAGS} -shared-intel")
endif ()
SaveCompilerFlags(C)
endif ()
