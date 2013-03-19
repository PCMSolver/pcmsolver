set(cavity_SOURCES 
    Atom.cpp
    Cavity.cpp
    GePolCavity.cpp
    Sphere.cpp
    SurfaceFunction.cpp
    WaveletCavity.cpp
)

add_library(cavity ${cavity_SOURCES})
add_dependencies(cavity external)

install(TARGETS cavity
 	ARCHIVE DESTINATION lib)
 
install(FILES Constants.h DESTINATION include)
install(FILES Cavity.h DESTINATION include)
install(FILES GePolCavity.h DESTINATION include)
install(FILES WaveletCavity.h DESTINATION include)
install(FILES Atom.h DESTINATION include)
install(FILES Sphere.h DESTINATION include)
