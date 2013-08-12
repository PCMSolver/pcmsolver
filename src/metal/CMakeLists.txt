set(metal_SOURCES
    gsfera.F
    plgndr.F
   )

add_library(metal ${metal_SOURCES})

# Export the sources up to src directory level
set(local_sources "")
foreach(source ${metal_SOURCES})
	get_filename_component(sourcepath ${source} ABSOLUTE)
  	set(local_sources ${sourcepath} ${local_sources})
endforeach()
set(SOURCES ${local_sources} ${SOURCES} PARENT_SCOPE)

install(TARGETS metal ARCHIVE DESTINATION lib)

