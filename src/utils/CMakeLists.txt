set(utils_SOURCES
	Input.cpp
	Solvent.cpp
	Atom.cpp
	Sphere.cpp	
    	SurfaceFunction.cpp
)

add_library(utils ${utils_SOURCES})
add_dependencies(utils external)

install(TARGETS utils
	ARCHIVE DESTINATION lib)

install(FILES PhysicalConstants.h DESTINATION include)
install(FILES Input.h DESTINATION include)
install(FILES Solvent.h DESTINATION include)
install(FILES Atom.h DESTINATION include)
install(FILES Sphere.h DESTINATION include)
install(FILES SurfaceFunction.h DESTINATION include)
install(FILES SurfaceMap.h DESTINATION include)
