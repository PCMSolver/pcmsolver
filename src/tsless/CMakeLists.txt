set(tsless_SOURCES
	tsless_cavity_interface.F
  )

if(ENABLE_TESTS)
	add_library(tsless ${tsless_SOURCES})
	install(TARGETS tsless ARCHIVE DESTINATION lib)
endif()

# Export the sources up to src directory level
set(local_sources "")
foreach(source ${tsless_SOURCES})
	get_filename_component(sourcepath ${source} ABSOLUTE)
  	set(local_sources ${sourcepath} ${local_sources})
endforeach()
set(SOURCES ${local_sources} ${SOURCES} PARENT_SCOPE)
