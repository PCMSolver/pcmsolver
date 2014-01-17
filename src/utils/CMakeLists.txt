set(utils_SOURCES
	Input.cpp
	Solvent.cpp
	Atom.cpp
	Sphere.cpp	
    	SurfaceFunction.cpp
	cnpy.cpp
)

if(ENABLE_TESTS)
	add_library(utils ${utils_SOURCES})
	add_dependencies(utils libgetkw)
	install(TARGETS utils ARCHIVE DESTINATION lib)
endif()

# Export the sources up to src directory level
set(local_sources "")
foreach(source ${utils_SOURCES})
	get_filename_component(sourcepath ${source} ABSOLUTE)
  	set(local_sources ${sourcepath} ${local_sources})
endforeach()
set(SOURCES ${local_sources} ${SOURCES} PARENT_SCOPE)

install(FILES PhysicalConstants.hpp DESTINATION include)
install(FILES Input.hpp DESTINATION include)
install(FILES Solvent.hpp DESTINATION include)
install(FILES Atom.hpp DESTINATION include)
install(FILES Sphere.hpp DESTINATION include)
install(FILES SurfaceFunction.hpp DESTINATION include)
install(FILES Citation.hpp DESTINATION include)
install(FILES CavityData.hpp DESTINATION include)
install(FILES GreenData.hpp DESTINATION include)
install(FILES SolverData.hpp DESTINATION include)
install(FILES cnpy.hpp DESTINATION include)
install(FILES MathUtils.hpp DESTINATION include)
