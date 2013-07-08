set(utils_SOURCES
	Input.cpp
	Solvent.cpp
	Atom.cpp
	Sphere.cpp	
    	SurfaceFunction.cpp
)

add_library(utils ${utils_SOURCES})

add_dependencies(utils libgetkw)

install(TARGETS utils ARCHIVE DESTINATION lib)

install(FILES PhysicalConstants.hpp DESTINATION include)
install(FILES Input.hpp DESTINATION include)
install(FILES Solvent.hpp DESTINATION include)
install(FILES Atom.hpp DESTINATION include)
install(FILES Sphere.hpp DESTINATION include)
install(FILES SurfaceFunction.hpp DESTINATION include)
