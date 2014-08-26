set(sources_list "")
# Get sources with Fortran extension
get_fortran_sources(sources_list)
# If you want to remove some sources specify them explictly here
list(REMOVE_ITEM sources_list "")

if(ENABLE_TESTS)
	add_library(metal STATIC ${sources_list})
	set(LIBS_TO_MERGE metal ${LIBS_TO_MERGE} PARENT_SCOPE)
	install(TARGETS metal ARCHIVE DESTINATION lib)
endif()

# Export the sources up to src directory level
set(local_sources "")
foreach(source ${sources_list})
	get_filename_component(sourcepath ${source} ABSOLUTE)
  	set(local_sources ${sourcepath} ${local_sources})
endforeach()
set(SOURCES ${local_sources} ${SOURCES} PARENT_SCOPE)

