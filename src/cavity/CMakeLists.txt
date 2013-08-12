set(cavity_SOURCES 
    GePolCavity.cpp
    WaveletCavity.cpp
    CavityFactory.cpp
)

add_library(cavity ${cavity_SOURCES})

# Export the sources up to src directory level
set(local_sources "")
foreach(source ${cavity_SOURCES})
	get_filename_component(sourcepath ${source} ABSOLUTE)
  	set(local_sources ${sourcepath} ${local_sources})
endforeach()
set(SOURCES ${local_sources} ${SOURCES} PARENT_SCOPE)

install(TARGETS cavity ARCHIVE DESTINATION lib)
 
install(FILES Cavity.hpp DESTINATION include)
install(FILES GePolCavity.hpp DESTINATION include)
install(FILES WaveletCavity.hpp DESTINATION include)
install(FILES CavityFactory.hpp DESTINATION include)
