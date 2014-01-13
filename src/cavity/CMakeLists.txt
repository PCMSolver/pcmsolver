set(cavity_SOURCES
    Cavity.cpp 
    GePolCavity.cpp 
    TsLessCavity.cpp
    WaveletCavity.cpp
    CavityFactory.cpp
)

if(ENABLE_TESTS)
	add_library(cavity ${cavity_SOURCES})
    if(DEVELOPMENT_CODE)
 	    add_dependencies(cavity pedra tsless wavcav)
    else()
 	    add_dependencies(cavity pedra tsless)
    endif()
	install(TARGETS cavity ARCHIVE DESTINATION lib)
endif()

# Export the sources up to src directory level
set(local_sources "")
foreach(source ${cavity_SOURCES})
	get_filename_component(sourcepath ${source} ABSOLUTE)
  	set(local_sources ${sourcepath} ${local_sources})
endforeach()
set(SOURCES ${local_sources} ${SOURCES} PARENT_SCOPE)
 
install(FILES Cavity.hpp DESTINATION include)
install(FILES GePolCavity.hpp DESTINATION include)
install(FILES TsLessCavity.hpp DESTINATION include)
install(FILES WaveletCavity.hpp DESTINATION include)
install(FILES CavityFactory.hpp DESTINATION include)
