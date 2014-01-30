set(pedra_SOURCES 
  pedra_utils.F90
  pedra_print.F90
  cavity_derivatives.F90
  pedra_cavity.F90
  pedra_cavity_interface.F90
  pedra_dblas.F90
  pedra_dlapack.F90
  makecav.F90
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
