set(pedra_SOURCES 
  pcm_utils.F
  pcm_print.F
  pcm_inertia.F90
  cavity_derivatives.F
  pedra_cavity.F
  pedra_cavity_interface.F
  makecav.F
  pcm_dblas1.F
  pcm_dlapack.F
  )

if(ENABLE_TESTS)
	add_library(pedra ${pedra_SOURCES})
	install(TARGETS pedra ARCHIVE DESTINATION lib)
endif()

# Export the sources up to src directory level
set(local_sources "")
foreach(source ${pedra_SOURCES})
	get_filename_component(sourcepath ${source} ABSOLUTE)
  	set(local_sources ${sourcepath} ${local_sources})
endforeach()
set(SOURCES ${local_sources} ${SOURCES} PARENT_SCOPE)
